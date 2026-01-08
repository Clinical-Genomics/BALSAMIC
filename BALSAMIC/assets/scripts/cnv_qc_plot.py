#!/usr/bin/env python

import math
from pathlib import Path

import click
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colormaps


# -----------------------------
# Helpers
# -----------------------------
def extract_info_float(info_str: str, key: str) -> float:
    """Extract float like 'MQ=60.00' from INFO; NaN if missing."""
    if pd.isna(info_str):
        return math.nan
    key_eq = key + "="
    for field in str(info_str).split(";"):
        if field.startswith(key_eq):
            try:
                return float(field[len(key_eq) :])
            except ValueError:
                return math.nan
    return math.nan


def parse_sample_fields(format_str: str, sample_str: str):
    """
    Parse AD and DP from FORMAT/TUMOR and compute BAF = alt / (ref + alt).

    No GT-based logic – we always try to compute BAF if AD is present.
    Returns (BAF, dummy_GT, DP_sample).
    """
    if pd.isna(format_str) or pd.isna(sample_str):
        return math.nan, None, math.nan

    f_keys = format_str.split(":")
    f_vals = sample_str.split(":")
    fmt = dict(zip(f_keys, f_vals))

    # DP from FORMAT
    dp_str = fmt.get("DP", None)
    try:
        dp_sample = int(dp_str) if dp_str not in (None, ".", "") else math.nan
    except ValueError:
        dp_sample = math.nan

    # AD: ref,alt(,alt2,...) – we use ref and first alt
    ad = fmt.get("AD", None)
    if ad is None or ad == ".":
        return math.nan, None, dp_sample

    try:
        counts = [int(x) for x in ad.split(",")]
        if len(counts) < 2:
            return math.nan, None, dp_sample
        ref_c = counts[0]
        alt_c = counts[1]  # if multi-allelic, this is the first alt
        total = ref_c + alt_c
        if total == 0:
            return math.nan, None, dp_sample
        baf = alt_c / total
        return baf, None, dp_sample
    except Exception:
        return math.nan, None, dp_sample


def flag_cnv_genes(df: pd.DataFrame) -> pd.DataFrame:
    """
    Subset PureCN genes to those that look like CNV genes.

    - drop 'Antitarget' and '-' pseudo-genes
    - keep genes where:
        C != 2 OR loh == TRUE OR type != NA
    """
    df = df.copy()
    df = df[~df["gene.symbol"].isin(["Antitarget", "-"])]

    cnv_mask = (
        (df["C"].notna()) & (df["C"] != 2)
    ) | (
        df["loh"].astype(str).str.upper() == "TRUE"
    ) | (
        df["type"].notna() & (df["type"].astype(str) != "NA")
    )

    return df[cnv_mask]


def normalize_chr_column(df: pd.DataFrame, prefer_col: str | None = None):
    """
    Detect and normalize a chromosome column to plain '1'..'22','X','Y'.
    Returns (df_with_chr, chr_col_name).
    """
    df = df.copy()

    candidate_cols = []
    if prefer_col and prefer_col in df.columns:
        candidate_cols.append(prefer_col)
    candidate_cols += [c for c in ["chr", "chromosome", "CHROM"] if c in df.columns]

    if not candidate_cols:
        raise ValueError("Could not find a chromosome column in dataframe.")

    col = candidate_cols[0]
    df[col] = df[col].astype(str)
    df[col] = df[col].str.replace("^chr", "", regex=True)

    return df, col


# -----------------------------
# Main plotting logic
# -----------------------------
def plot_chromosomes(
    pon_path: Path,
    cnr_path: Path,
    vcf_path: Path,
    genes_path: Path,
    segs_path: Path,
    outdir: Path,
    include_y: bool = False,
    weight_thresh: float = 0.5,
    pct_spread: float = 0.90,
    window: int = 21,
    anti_factor: float = 0.2,
    base_label_offset: float = 0.2,
    dp_min: int = 15,
    mq_min: float = 35.0,
    qd_min: float = 30.0,
    sor_max: float = 2.0,
):
    outdir.mkdir(parents=True, exist_ok=True)

    # ----------------------------
    # Chromosome order
    # ----------------------------
    chr_order = [str(i) for i in range(1, 23)] + ["X"]
    if include_y:
        chr_order.append("Y")

    # ----------------------------
    # 1. Load PON and CNR
    # ----------------------------
    pon = pd.read_csv(pon_path, sep="\t")
    cnr = pd.read_csv(cnr_path, sep="\t")

    pon["chromosome"] = pon["chromosome"].astype(str)
    cnr["chromosome"] = cnr["chromosome"].astype(str)

    pon = pon.dropna(subset=["spread"])

    pon = pon[pon["chromosome"].isin(chr_order)]
    cnr = cnr[cnr["chromosome"].isin(chr_order)]

    # 2. Merge PON and CNR on bins
    merged = pd.merge(
        cnr,
        pon[["chromosome", "start", "end", "spread"]],
        on=["chromosome", "start", "end"],
        how="inner",
        suffixes=("_cnr", "_pon"),
    )

    # ----------------------------
    # 3. Load genes & segments from PureCN
    # ----------------------------
    genes = pd.read_csv(genes_path)
    genes, genes_chr_col = normalize_chr_column(genes)
    genes = genes[genes[genes_chr_col].isin(chr_order)]

    genes_cnv = flag_cnv_genes(genes)
    # require at least 4 targets for label coloring (tune if you like)
    genes_cnv = genes_cnv[genes_cnv["number.targets"] > 4]

    segs = pd.read_csv(segs_path)
    segs, segs_chr_col = normalize_chr_column(segs)
    segs = segs[segs[segs_chr_col].isin(chr_order)]

    # ----------------------------
    # 4. Load VCF and compute BAF, DP, MQ, QD, SOR
    # ----------------------------
    vcf_cols = [
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO",
        "FORMAT",
        "TUMOR",
    ]
    vcf = pd.read_csv(vcf_path, sep="\t", comment="#", header=None, names=vcf_cols)

    vcf["CHROM"] = vcf["CHROM"].astype(str)
    vcf["CHROM"] = vcf["CHROM"].str.replace("^chr", "", regex=True)
    vcf = vcf[vcf["CHROM"].isin(chr_order)]

    bafs, gts, dps = zip(
        *[
            parse_sample_fields(fmt, sample)
            for fmt, sample in zip(vcf["FORMAT"], vcf["TUMOR"])
        ]
    )

    vcf["BAF"] = bafs
    vcf["DP_sample"] = dps

    vcf["MQ"] = vcf["INFO"].apply(lambda s: extract_info_float(s, "MQ"))
    vcf["QD"] = vcf["INFO"].apply(lambda s: extract_info_float(s, "QD"))
    vcf["SOR"] = vcf["INFO"].apply(lambda s: extract_info_float(s, "SOR"))

    # Quality filters (no BAF-based filtering)
    vcf = vcf[
        (vcf["DP_sample"] >= dp_min)
        & (vcf["MQ"] > mq_min)
        & (vcf["QD"] > qd_min)
        & (vcf["SOR"] < sor_max)
    ]

    # ----------------------------
    # 5. QC thresholds & smoothing setup
    # ----------------------------
    spread_thresh = merged["spread"].quantile(pct_spread)

    y_max = max(merged["log2"].abs().max(), merged["spread"].max())
    y_pad = 0.1
    y_lim = (-y_max - y_pad, y_max + y_pad)

    # ----------------------------
    # 6. Per-chromosome plots
    # ----------------------------
    for chr_name in chr_order:
        sub = merged[merged["chromosome"] == chr_name].copy()
        if sub.empty:
            continue

        # Sort by genomic position
        sub = sub.sort_values("start")

        # Classify target vs antitarget
        sub["type"] = np.where(sub["gene"] == "Antitarget", "Antitarget", "Target")

        # Plotting filters
        sub = sub[(sub["weight"] > weight_thresh) & (sub["spread"] <= spread_thresh)]
        if sub.empty:
            continue

        # Re-sort (after filtering) and assign variable-width x-coordinate
        sub = sub.sort_values("start")
        sub["bin_width"] = np.where(sub["type"] == "Antitarget", anti_factor, 1.0)
        sub["x_coord"] = sub["bin_width"].cumsum() - sub["bin_width"] / 2

        # Smoothing (rolling on row order; x only used for plotting)
        sub = sub.sort_values("x_coord")
        sub["log2_smooth"] = sub["log2"].rolling(
            window=window, center=True
        ).median()
        sub["spread_smooth"] = sub["spread"].rolling(
            window=window, center=True
        ).median()

        # Helper: map genomic POS → x_coord (nearest bin)
        bin_starts = sub["start"].values
        x_coords = sub["x_coord"].values

        def pos_to_xcoord(pos):
            idx = np.searchsorted(bin_starts, pos, side="right") - 1
            if idx < 0:
                idx = 0
            if idx >= len(x_coords):
                idx = len(x_coords) - 1
            return x_coords[idx]

        # BAF for this chromosome
        baf_chr = vcf[vcf["CHROM"] == chr_name].copy()
        baf_chr = baf_chr.sort_values("POS")

        if not baf_chr.empty:
            baf_chr["x_coord"] = baf_chr["POS"].apply(pos_to_xcoord)

        # CNV genes for this chromosome
        genes_chr = genes_cnv[genes_cnv[genes_chr_col] == chr_name].copy()
        genes_chr = genes_chr.sort_values("start")

        gene_to_color: dict[str, tuple] = {}
        gene_names = np.array([])

        if not genes_chr.empty:
            genes_chr["x_coord"] = genes_chr["start"].apply(pos_to_xcoord)

            gene_names = genes_chr["gene.symbol"].unique()
            n_genes = len(gene_names)
            cmap = colormaps.get_cmap("tab20").resampled(n_genes)
            for i, g in enumerate(gene_names):
                gene_to_color[g] = cmap(i)

        # Segments for this chromosome
        segs_chr = segs[segs[segs_chr_col] == chr_name].copy()
        segs_chr = segs_chr.sort_values("start")

        if not segs_chr.empty:
            segs_chr["x_start"] = segs_chr["start"].apply(pos_to_xcoord)
            segs_chr["x_end"] = segs_chr["end"].apply(pos_to_xcoord)

        # ---------------- Figure + axes ----------------
        fig, (ax1, ax2) = plt.subplots(
            2,
            1,
            figsize=(14, 6),
            sharex=True,
            gridspec_kw={"height_ratios": [2, 1]},
        )

        x = sub["x_coord"]

        # Background split: non-CNV targets & antitargets
        if len(gene_names) > 0:
            non_cnv = sub[~sub["gene"].isin(gene_names)]
        else:
            non_cnv = sub

        bg_targets = non_cnv[non_cnv["type"] == "Target"]
        bg_antis = non_cnv[non_cnv["type"] == "Antitarget"]

        # Non-CNV antitargets: faint grey
        if not bg_antis.empty:
            ax1.scatter(
                bg_antis["x_coord"],
                bg_antis["log2"],
                s=3,
                alpha=0.10,
                color="lightgrey",
                label="Antitarget bins",
            )

        # Non-CNV targets: visible neutral
        if not bg_targets.empty:
            ax1.scatter(
                bg_targets["x_coord"],
                bg_targets["log2"],
                s=4,
                alpha=0.4,
                color="tab:blue",
                label="Target bins (no CNV gene)",
            )

        # CNV genes: colored points
        if len(gene_names) > 0:
            for gene in gene_names:
                gsub = sub[sub["gene"] == gene]
                if gsub.empty:
                    continue
                color = gene_to_color[gene]
                ax1.scatter(
                    gsub["x_coord"],
                    gsub["log2"],
                    s=8,
                    alpha=0.9,
                    color=color,
                )

        # Filled spread band (PON noise envelope)
        ax1.fill_between(
            x,
            -sub["spread_smooth"],
            sub["spread_smooth"],
            alpha=0.3,
            step="mid",
            label="PON noise band",
        )

        # Smoothed log2
        ax1.plot(
            x,
            sub["log2_smooth"],
            linewidth=1.5,
            alpha=0.9,
            color="tab:green",
            label=f"log2 (median {window} bins)",
        )

        # Segment bars
        if not segs_chr.empty:
            for _, srow in segs_chr.iterrows():
                xs = srow["x_start"]
                xe = srow["x_end"]
                y_seg = srow["seg.mean"]

                C_val = srow.get("C", math.nan)
                if pd.isna(C_val):
                    seg_color = "black"
                elif C_val > 2:
                    seg_color = "red"
                elif C_val < 2:
                    seg_color = "royalblue"
                else:
                    seg_color = "black"

                ax1.hlines(
                    y_seg,
                    xs,
                    xe,
                    colors=seg_color,
                    linewidth=1.0,
                    alpha=0.8,
                )

        ax1.axhline(0, color="black", linewidth=0.8)
        ax1.set_ylim(*y_lim)
        ax1.set_ylabel("log2 / spread")
        ax1.set_title(
            f"Chr {chr_name} – log2 vs PON spread + BAF + CNV genes + segments\n"
            "(targets expanded, antitargets compressed)"
        )
        ax1.legend(loc="upper right", fontsize=8)

        # Gene labels (directional offset)
        if not genes_chr.empty and len(gene_names) > 0:
            seen_positions = set()

            for _, row in genes_chr.iterrows():
                gene = row["gene.symbol"]
                gx = row["x_coord"]
                key = (gene, gx)
                if key in seen_positions:
                    continue
                seen_positions.add(key)

                color = gene_to_color.get(gene, "black")

                idx_near = np.argmin(np.abs(sub["x_coord"] - gx))
                y_val = sub["log2_smooth"].iloc[idx_near]
                if pd.isna(y_val):
                    y_val = 0.0

                # gains up, losses down
                offset = base_label_offset if y_val >= 0 else -base_label_offset
                y_label = y_val + offset

                ax1.axvline(gx, color=color, linewidth=0.5, alpha=0.3)
                ax1.text(
                    gx,
                    y_label,
                    gene,
                    rotation=90,
                    fontsize=8,
                    ha="center",
                    va="bottom" if offset > 0 else "top",
                    color=color,
                )

        # BAF panel
        if not baf_chr.empty:
            ax2.scatter(
                baf_chr["x_coord"],
                baf_chr["BAF"],
                s=6,
                alpha=0.5,
            )
        ax2.axhline(0.5, color="gray", linewidth=0.8, linestyle="--")
        for frac in [1 / 3, 2 / 3]:
            ax2.axhline(frac, color="lightgray", linewidth=0.6, linestyle=":")
        ax2.set_ylim(0, 1)
        ax2.set_xlabel(
            "Pseudo-position (targets expanded, antitargets compressed)"
        )
        ax2.set_ylabel("BAF")

        plt.tight_layout()
        out_path = outdir / f"qc_chr{chr_name}_segments.png"
        plt.savefig(out_path, dpi=150)
        plt.close(fig)


# -----------------------------
# Click CLI
# -----------------------------
@click.command()
@click.option("--pon", "pon_path", required=True, type=click.Path(exists=True), help="CNVkit PON .cnn file")
@click.option("--cnr", "cnr_path", required=True, type=click.Path(exists=True), help="Tumor CNVkit .cnr file")
@click.option("--vcf", "vcf_path", required=True, type=click.Path(exists=True), help="VCF with germline SNPs")
@click.option("--genes", "genes_path", required=True, type=click.Path(exists=True), help="PureCN *_genes.csv file")
@click.option("--segs", "segs_path", required=True, type=click.Path(exists=True), help="PureCN *_loh.csv or segments CSV")
@click.option("--outdir", required=True, type=click.Path(), help="Output directory for PNGs")
@click.option("--include-y/--no-include-y", default=False, show_default=True, help="Include chromosome Y")
@click.option("--weight-thresh", default=0.5, show_default=True, help="Minimum bin weight to keep")
@click.option("--spread-quantile", default=0.90, show_default=True, help="Spread quantile threshold for noisy bins")
@click.option("--window", default=21, show_default=True, help="Rolling window size for smoothing")
@click.option("--anti-factor", default=0.2, show_default=True, help="Relative x-width of antitarget bins")
@click.option("--label-offset", default=0.2, show_default=True, help="Vertical offset for gene labels (log2 units)")
@click.option("--dp-min", default=15, show_default=True, help="Minimum DP for SNPs used in BAF")
@click.option("--mq-min", default=35.0, show_default=True, help="Minimum MQ for SNPs used in BAF")
@click.option("--qd-min", default=30.0, show_default=True, help="Minimum QD for SNPs used in BAF")
@click.option("--sor-max", default=2.0, show_default=True, help="Maximum SOR for SNPs used in BAF")
def cli(
    pon_path,
    cnr_path,
    vcf_path,
    genes_path,
    segs_path,
    outdir,
    include_y,
    weight_thresh,
    spread_quantile,
    window,
    anti_factor,
    label_offset,
    dp_min,
    mq_min,
    qd_min,
    sor_max,
):
    """Generate per-chromosome QC plots: log2 + PON spread + BAF + PureCN genes/segments."""
    plot_chromosomes(
        pon_path=Path(pon_path),
        cnr_path=Path(cnr_path),
        vcf_path=Path(vcf_path),
        genes_path=Path(genes_path),
        segs_path=Path(segs_path),
        outdir=Path(outdir),
        include_y=include_y,
        weight_thresh=weight_thresh,
        pct_spread=spread_quantile,
        window=window,
        anti_factor=anti_factor,
        base_label_offset=label_offset,
        dp_min=dp_min,
        mq_min=mq_min,
        qd_min=qd_min,
        sor_max=sor_max,
    )


if __name__ == "__main__":
    cli()
