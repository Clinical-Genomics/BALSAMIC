import math
from pathlib import Path
from typing import Dict, Tuple, List
from pandas.errors import EmptyDataError
import fitz
import matplotlib.pyplot as plt
from matplotlib import colormaps
import numpy as np
import pandas as pd


# =============================================================================
# Generic helpers
# =============================================================================


def pdf_first_page_to_png(pdf_path: str, png_path: str, dpi: int = 300) -> None:
    doc = fitz.open(pdf_path)
    page = doc[0]
    page.get_pixmap(dpi=dpi).save(png_path)
    doc.close()


def load_cancer_gene_set(
    path: str | Path, min_occurrence: int = 1, only_annotated: bool = True
) -> set[str]:
    df = pd.read_csv(path, sep="\t", dtype=str)
    # Clean up column names
    df.columns = df.columns.str.strip()

    # Standardize column names we care about
    # (Change these strings if your headers differ slightly.)
    symbol_col = "Hugo Symbol"
    occ_col = "# of occurrence within resources (Column J-P)"
    onco_col = "OncoKB Annotated"
    type_col = "Gene Type"

    # Convert occurrence to numeric
    df[occ_col] = pd.to_numeric(df[occ_col], errors="coerce").fillna(0).astype(int)

    # Basic filters: cancer-relevant genes
    mask = df[occ_col] >= min_occurrence

    if only_annotated and onco_col in df.columns:
        mask &= df[onco_col].str.upper().eq("YES")

    # Optional: keep classic cancer gene types only
    if type_col in df.columns:
        mask &= df[type_col].isin(["ONCOGENE", "TSG"])

    selected = df.loc[mask, symbol_col].astype(str).str.strip()
    return set(selected)


def normalize_chr_column(
    df: pd.DataFrame, prefer_col: str | None = None
) -> tuple[pd.DataFrame, str]:
    """
    Detect and normalize a chromosome column to plain '1'..'22','X','Y'.
    Returns (df_with_chr, chr_col_name).
    """
    df = df.copy()

    candidate_cols: list[str] = []
    if prefer_col and prefer_col in df.columns:
        candidate_cols.append(prefer_col)
    candidate_cols += [c for c in ["chr", "chromosome", "CHROM"] if c in df.columns]

    if not candidate_cols:
        raise ValueError("Could not find a chromosome column in dataframe.")

    col = candidate_cols[0]
    df[col] = df[col].astype(str)
    df[col] = df[col].str.replace("^chr", "", regex=True)

    return df, col


def parse_sample_fields(format_str: str, sample_str: str) -> tuple[float, None, float]:
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


def safe_read_csv(path: str | Path, **kwargs) -> pd.DataFrame:
    """
    Wrapper around pd.read_csv that returns an empty DataFrame
    if file is empty or missing.
    """
    try:
        return pd.read_csv(path, **kwargs)
    except (EmptyDataError, FileNotFoundError):
        return pd.DataFrame()


# =============================================================================
# Cytoband
# =============================================================================


def load_cytobands(path: str | Path) -> pd.DataFrame:
    cols = ["chrom", "chromStart", "chromEnd", "name", "gieStain"]
    cyto = pd.read_csv(path, sep="\t", header=None, names=cols)

    # normalize
    cyto["chrom"] = cyto["chrom"].astype(str).str.replace("^chr", "", regex=True)
    cyto["start_int"] = cyto["chromStart"].astype(int)
    cyto["end_int"] = cyto["chromEnd"].astype(int)

    return cyto


def annotate_genes_with_cytoband(
    df_genes: pd.DataFrame, cyto: pd.DataFrame
) -> pd.DataFrame:
    """
    Add a 'cytoband' column to df_genes based on the segment coordinates
    (seg_start, seg_end). If seg_* are NaN, falls back to gene start/end.

    If a segment overlaps multiple bands, returns a range like '7q21.1-q22.3'.
    If exactly one band, returns e.g. '7q22.1'.
    """
    df = df_genes.copy()

    if "chr" not in df.columns:
        raise ValueError("df_genes must contain a 'chr' column.")

    # normalize chr to match cytoband
    df["chr"] = df["chr"].astype(str).str.replace("^chr", "", regex=True)

    # Prepare output column
    df["cytoband"] = pd.Series(pd.NA, index=df.index, dtype="string")

    # ensure seg_start / seg_end exist; if not, create as NaN
    if "seg_start" not in df.columns:
        df["seg_start"] = np.nan
    if "seg_end" not in df.columns:
        df["seg_end"] = np.nan

    for chrom, genes_chr in df.groupby("chr"):
        cyto_chr = cyto[cyto["chrom"] == chrom]
        if cyto_chr.empty:
            continue

        cyto_chr = cyto_chr.sort_values("start_int")
        band_starts = cyto_chr["start_int"].to_numpy()
        band_ends = cyto_chr["end_int"].to_numpy()
        band_names = cyto_chr["name"].to_numpy()

        for idx, row in genes_chr.iterrows():
            # prefer segment coords; fallback to gene coords
            if not pd.isna(row["seg_start"]) and not pd.isna(row["seg_end"]):
                s_start = int(row["seg_start"])
                s_end = int(row["seg_end"])
            else:
                if "start" not in df.columns or "end" not in df.columns:
                    continue
                s_start = int(row["start"])
                s_end = int(row["end"])

            mask = (band_starts <= s_end) & (band_ends >= s_start)
            if not mask.any():
                continue

            names = band_names[mask]

            if len(names) == 1:
                label = f"{chrom}{names[0]}"
            else:
                first = names[0]
                last = names[-1]
                label = f"{chrom}{first}-{chrom}{last}"

            df.at[idx, "cytoband"] = label

    return df


# =============================================================================
# PureCN CNV genes (for plotting only)
# =============================================================================


def flag_cnv_genes(df: pd.DataFrame, min_targets_cnvgene: int = 4) -> pd.DataFrame:
    """
    Subset PureCN genes to those that look like CNV genes.

      - drop 'Antitarget' and '-' pseudo-genes
      - keep genes where:
          C != 2 OR loh == TRUE OR type != NA
      - require at least min_targets_cnvgene targets
    """
    df = df.copy()
    df = df[~df["gene.symbol"].isin(["Antitarget", "-"])]

    cnv_mask = (
        ((df["C"].notna()) & (df["C"] != 2))
        | (df["loh"].astype(str).str.upper() == "TRUE")
        | (df["type"].notna() & (df["type"].astype(str) != "NA"))
    )

    df = df[cnv_mask]

    if "number.targets" in df.columns and min_targets_cnvgene is not None:
        df = df[df["number.targets"] > min_targets_cnvgene]

    return df


# =============================================================================
# CNR + PON merge
# =============================================================================


def merge_cnr_pon(
    df_cnr: pd.DataFrame,
    df_pon: pd.DataFrame,
    chr_col: str = "chromosome",
    spread_col: str = "spread",
) -> pd.DataFrame:
    """
    Merge CNVkit tumor CNR with PON on (chr, start, end) and drop rows
    without PON spread.
    """
    cnr = df_cnr.copy()
    pon = df_pon.copy()

    cnr[chr_col] = cnr[chr_col].astype(str).str.replace("^chr", "", regex=True)
    pon[chr_col] = pon[chr_col].astype(str).str.replace("^chr", "", regex=True)

    if spread_col not in pon.columns:
        raise ValueError(f"PON is missing required column '{spread_col}'")

    pon = pon.dropna(subset=[spread_col])

    merged = pd.merge(
        cnr,
        pon[[chr_col, "start", "end", spread_col]],
        on=[chr_col, "start", "end"],
        how="inner",
        suffixes=("_cnr", "_pon"),
    )

    # If PON spread came in as spread_pon, normalize back to 'spread'
    if spread_col + "_pon" in merged.columns:
        merged[spread_col] = merged[spread_col + "_pon"]

    return merged


# =============================================================================
# VCF → BAF
# =============================================================================


def load_vcf_with_baf(vcf_path: str | Path, chr_order: list[str]) -> pd.DataFrame:
    """
    Load a (possibly unfiltered) VCF with a single tumor sample and compute:

      - BAF (alt / (ref+alt))
      - DP_sample
    """
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

    bafs, _, dps = zip(
        *[
            parse_sample_fields(fmt, sample)
            for fmt, sample in zip(vcf["FORMAT"], vcf["TUMOR"])
        ]
    )

    vcf["BAF"] = bafs
    vcf["DP_sample"] = dps

    return vcf


# =============================================================================
# Build per-gene table from CNR
# =============================================================================


def build_genes_from_cnr(cnr_path: str | Path) -> pd.DataFrame:
    """
    Build a per-gene table from CNVkit .cnr:

      - one row per (chromosome, gene.symbol)
      - gene span from min(start)/max(end)
      - number.targets = number of bins for that gene
      - gene.mean / gene.min / gene.max from CNR log2

    Pseudo-genes 'Antitarget' and '-' are dropped.
    """
    cnr = safe_read_csv(cnr_path, sep="\t").copy()
    if cnr.empty:
        return cnr

    # normalize chromosome
    chr_col = "chromosome" if "chromosome" in cnr.columns else "chr"
    cnr[chr_col] = cnr[chr_col].astype(str).str.replace("^chr", "", regex=True)

    gene_col = "gene"
    if gene_col not in cnr.columns:
        raise ValueError("CNR file must contain a 'gene' column.")

    # drop pseudo-genes
    cnr = cnr[cnr[gene_col].notna()]
    cnr = cnr[~cnr[gene_col].isin(["Antitarget", "-"])]

    # expand multi-gene bins (A,B,C → 3 rows)
    cnr["gene.symbol"] = (
        cnr[gene_col]
        .astype(str)
        .str.split(",")
        .apply(lambda lst: [g.strip() for g in lst if g.strip()])
    )
    cnr = cnr.explode("gene.symbol")

    # basic per-gene aggregation
    agg = cnr.groupby([chr_col, "gene.symbol"], as_index=False).agg(
        start=("start", "min"),
        end=("end", "max"),
        number_targets=("log2", "count"),
        gene_mean=("log2", "mean"),
        gene_min=("log2", "min"),
        gene_max=("log2", "max"),
        weight_mean=("weight", "mean"),
    )

    # match existing naming in your code
    agg = agg.rename(
        columns={
            chr_col: "chr",
            "number_targets": "number.targets",
            "gene_mean": "gene.mean",
            "gene_min": "gene.min",
            "gene_max": "gene.max",
            "weight_mean": "weight.mean",
        }
    )

    # round gene stats a bit
    for col in ["gene.mean", "gene.min", "gene.max", "weight.mean"]:
        agg[col] = agg[col].round(3)

    return agg


# =============================================================================
# Per-chromosome plots (PON spread + log2 + segments + BAF + CNV genes)
# =============================================================================


def plot_chromosomes(
    pon_path: Path | None,
    cnr_path: Path,
    cns_path: Path,
    vcf_path: Path,
    genes_path: Path,
    segs_path: Path,
    outdir: Path,
    case_id: str,
    include_y: bool = False,
    weight_thresh: float = 0.1,
    pct_spread: float = 0.95,
    window: int = 5,
    anti_factor: float = 0.15,
    base_label_offset: float = 1.5,
    min_targets_cnvgene: int = 2,
    cancer_genes: set[str] | None = None,
    neutral_target_factor: float = 0.3,
    is_exome: bool = False,
) -> None:
    """
    Main plotting routine creating one PNG per chromosome with:

      - PON spread band (if PON available)
      - smoothed log2 coverage
      - CNV segments (PureCN LOH regions)
      - CNV genes colored & labeled (PureCN genes)
      - BAF panel underneath

    For large panels / exomes:
      - If `cancer_genes` is provided, only CNV genes whose symbol is in this
        set are highlighted & labeled.
      - Non-CNV targets are horizontally compressed by `neutral_target_factor`
        so that cancer-relevant CNV genes stand out.

    Y-axis is scaled *per chromosome* (log2, spread, segments) instead of
    globally across all chromosomes.
    """
    outdir.mkdir(parents=True, exist_ok=True)

    # Chromosome order
    chr_order = [str(i) for i in range(1, 23)] + ["X"]
    if include_y:
        chr_order.append("Y")

    # 1. Load CNR (always) and PON (optional) and merge
    cnr = safe_read_csv(cnr_path, sep="\t")
    if cnr.empty:
        return

    use_pon = False
    merged = None

    if pon_path is not None and Path(pon_path).is_file():
        pon = safe_read_csv(pon_path, sep="\t")
        if not pon.empty and "spread" in pon.columns:
            try:
                merged = merge_cnr_pon(
                    cnr,
                    pon,
                    chr_col="chromosome",
                    spread_col="spread",
                )
                if not merged.empty:
                    use_pon = True
            except Exception:
                use_pon = False

    if not use_pon:
        # Fall back: use cnr only, with a dummy spread column so code still works
        merged = cnr.copy()
        merged["chromosome"] = (
            merged["chromosome"].astype(str).str.replace("^chr", "", regex=True)
        )
        if "spread" not in merged.columns:
            merged["spread"] = np.nan

    merged = merged[merged["chromosome"].isin(chr_order)]

    # 2. Load genes & segments from PureCN
    genes = safe_read_csv(genes_path)
    genes_cnv = pd.DataFrame()
    genes_chr_col = "chr"
    if not genes.empty:
        genes, genes_chr_col = normalize_chr_column(genes)
        genes = genes[genes[genes_chr_col].isin(chr_order)]

        # --- Pass 1: standard CNV gene calling (e.g. needs >= 2 targets) ---
        genes_cnv = flag_cnv_genes(
            genes,
            min_targets_cnvgene=min_targets_cnvgene,
        )

        # --- Pass 2: relaxed calling for cancer genes (>= 1 target) ---
        if cancer_genes:
            # Work on the subset of genes that are cancer-relevant
            genes_cancer = genes[genes["gene.symbol"].isin(cancer_genes)].copy()
            if not genes_cancer.empty:
                genes_cnv_cancer_relaxed = flag_cnv_genes(
                    genes_cancer,
                    min_targets_cnvgene=1,  # <- always allow 1 target for cancer genes
                )

                # Union of (normal CNV genes) ∪ (relaxed cancer CNV genes)
                genes_cnv = pd.concat(
                    [genes_cnv, genes_cnv_cancer_relaxed], ignore_index=True
                ).drop_duplicates(subset=[genes_chr_col, "gene.symbol"])

            # For exome runs, we still want to *only* keep cancer genes
            if is_exome:
                genes_cnv = genes_cnv[genes_cnv["gene.symbol"].isin(cancer_genes)]

    segs = safe_read_csv(segs_path)
    segs_chr_col = "chr"
    if not segs.empty:
        segs, segs_chr_col = normalize_chr_column(segs)
        segs = segs[segs[segs_chr_col].isin(chr_order)]

    # 3. Load VCF and compute BAF
    vcf = load_vcf_with_baf(
        vcf_path=vcf_path,
        chr_order=chr_order,
    )

    # 4. QC thresholds (global spread threshold, per-chr plotting)
    if use_pon:
        spread_thresh = merged["spread"].quantile(pct_spread)
    else:
        spread_thresh = np.inf

    # 5. Per-chromosome plots
    for chr_name in chr_order:
        sub = merged[merged["chromosome"] == chr_name].copy()
        if sub.empty:
            continue

        # Sort by genomic position
        sub = sub.sort_values("start")

        # Classify target vs antitarget
        sub["type"] = np.where(sub["gene"] == "Antitarget", "Antitarget", "Target")

        # Plotting filters
        if use_pon:
            sub = sub[
                (sub["weight"] > weight_thresh) & (sub["spread"] <= spread_thresh)
            ]
        else:
            sub = sub[sub["weight"] > weight_thresh]
        if sub.empty:
            continue

        # Determine which CNV genes exist on this chromosome,
        # optionally restricted by `cancer_genes`.
        chr_genes = pd.DataFrame()
        gene_to_color: dict[str, tuple] = {}
        gene_names = np.array([])

        if not genes_cnv.empty:
            chr_genes = genes_cnv[genes_cnv[genes_chr_col] == chr_name].copy()
            chr_genes = chr_genes.sort_values("start")
            if not chr_genes.empty:
                gene_names = chr_genes["gene.symbol"].unique()

        # Mark which bins belong to CNV genes
        sub["is_cnv_gene"] = sub["gene"].isin(gene_names)

        # variable-width x-coordinate:
        # - antitargets compressed by anti_factor
        # - CNV-gene targets use width 1.0
        # - non-CNV targets compressed by neutral_target_factor
        def _bin_width(row: pd.Series) -> float:
            if row["type"] == "Antitarget":
                return anti_factor
            if row["is_cnv_gene"]:
                return 1.0
            return neutral_target_factor

        sub["bin_width"] = sub.apply(_bin_width, axis=1)
        sub["x_coord"] = sub["bin_width"].cumsum() - sub["bin_width"] / 2

        # smoothing
        sub = sub.sort_values("x_coord")
        sub["log2_smooth"] = sub["log2"].rolling(window=window, center=True).median()
        if use_pon:
            sub["spread_smooth"] = (
                sub["spread"].rolling(window=window, center=True).median()
            )
        else:
            sub["spread_smooth"] = np.nan

        # helper: map genomic POS → x_coord (nearest bin)
        bin_starts = sub["start"].values
        x_coords = sub["x_coord"].values

        def pos_to_xcoord(pos: int) -> float:
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

        # Build gene color map and place gene x-coordinates
        if not chr_genes.empty and len(gene_names) > 0:
            chr_genes["x_coord"] = chr_genes["start"].apply(pos_to_xcoord)
            n_genes = len(gene_names)
            cmap = colormaps.get_cmap("tab20").resampled(n_genes)
            for i, g in enumerate(gene_names):
                gene_to_color[g] = cmap(i)

        genes_chr = chr_genes  # keep old name for the labeling code

        # Segments for this chromosome (PureCN regions)
        segs_chr = pd.DataFrame()
        if not segs.empty:
            segs_chr = segs[segs[segs_chr_col] == chr_name].copy()
            segs_chr = segs_chr.sort_values("start")
            if not segs_chr.empty:
                segs_chr["x_start"] = segs_chr["start"].apply(pos_to_xcoord)
                segs_chr["x_end"] = segs_chr["end"].apply(pos_to_xcoord)

        # -------- Per-chromosome y-scaling --------
        # log2 range on this chromosome
        log2_max = sub["log2"].abs().max()

        # spread range on this chromosome (if PON)
        if use_pon and sub["spread_smooth"].notna().any():
            spread_max = sub["spread_smooth"].abs().max()
        elif use_pon:
            spread_max = sub["spread"].abs().max()
        else:
            spread_max = 0.0

        # segment range on this chromosome
        if not segs_chr.empty and "seg.mean" in segs_chr.columns:
            seg_max = segs_chr["seg.mean"].abs().max()
        else:
            seg_max = 0.0

        # Protect against NaNs and too-small ranges
        y_max = max(0.5, float(log2_max), float(spread_max), float(seg_max))
        y_pad = 0.1
        y_lim_chr = (-y_max - y_pad, y_max + y_pad)

        # ---------- Figure + axes ----------
        fig, (ax1, ax2) = plt.subplots(
            2,
            1,
            figsize=(14, 6),
            sharex=True,
            gridspec_kw={"height_ratios": [2, 1]},
        )

        x = sub["x_coord"]

        # background split: non-CNV targets & antitargets
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
                alpha=0.38,
                color="lightgrey",
                label="Antitarget bins",
            )

        # Non-CNV targets: visible but compressed neutrals
        if not bg_targets.empty:
            ax1.scatter(
                bg_targets["x_coord"],
                bg_targets["log2"],
                s=4,
                alpha=0.5,
                color="tab:blue",
                label="Target bins (no highlighted CNV gene)",
            )

        # CNV genes: colored points
        if len(gene_names) > 0:
            for gene in gene_names:
                gsub = sub[sub["gene"] == gene]
                if gsub.empty:
                    continue
                color = gene_to_color.get(gene, "black")
                ax1.scatter(
                    gsub["x_coord"],
                    gsub["log2"],
                    s=8,
                    alpha=0.9,
                    color=color,
                )

        # PON spread band
        if use_pon and sub["spread_smooth"].notna().any():
            ax1.fill_between(
                x,
                -sub["spread_smooth"],
                sub["spread_smooth"],
                alpha=0.4,
                step="mid",
                label="PON noise band",
            )

        # smoothed log2
        ax1.plot(
            x,
            sub["log2_smooth"],
            linewidth=1.5,
            alpha=0.9,
            color="tab:green",
            label=f"log2 (median {window} bins)",
        )

        # segment bars
        if not segs_chr.empty and "seg.mean" in segs_chr.columns:
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
        ax1.set_ylim(*y_lim_chr)
        ax1.set_ylabel("log2 / spread")
        title_suffix = "log2 vs PON spread" if use_pon else "log2 (no PON available)"
        ax1.set_title(f"Chr {chr_name} – {title_suffix}: {case_id}\n")
        ax1.legend(loc="upper right", fontsize=8)

        # Gene labels (directional offset)
        if not genes_chr.empty and len(gene_names) > 0:
            seen_positions: set[tuple[str, float]] = set()

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
                offset = base_label_offset if y_val < 0 else -base_label_offset
                y_label = y_val + offset

                ax1.axvline(gx, color=color, linewidth=0.5, alpha=0.3)
                ax1.text(
                    gx,
                    y_label,
                    gene,
                    rotation=90,
                    fontsize=10,
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
        ax2.set_xlabel("Pseudo-position (CNV genes expanded, other bins compressed)")
        ax2.set_ylabel("BAF")

        plt.tight_layout()
        out_path = outdir / f"cnv_chr{chr_name}_segments.png"
        plt.savefig(out_path, dpi=150)
        plt.close(fig)


# =============================================================================
# Segment → gene table annotation (size, MAF, seg.mean, num.snps, C/M/loh/type)
# =============================================================================


def annotate_genes_with_segments(
    df_genes: pd.DataFrame, df_regions: pd.DataFrame
) -> pd.DataFrame:
    """
    For each gene row in df_genes, find the LOH region in df_regions on the same
    chromosome with the largest overlap, and annotate:

        seg_start, seg_end, seg_size_bp,
        maf.expected, maf.observed (if present in df_regions),
        seg.mean, num.snps, C, M, M.flagged (if present)

    Matching is purely by chromosome and genomic intervals (no seg.id needed).
    """
    df_genes = df_genes.copy()
    df_regions = df_regions.copy()
    if df_regions.empty:
        # ensure seg_* columns exist
        for col in ["seg_start", "seg_end", "seg_size_bp"]:
            if col not in df_genes.columns:
                df_genes[col] = np.nan
        return df_genes

    # normalize chr / positions
    df_genes["chr"] = df_genes["chr"].astype(str).str.replace("^chr", "", regex=True)
    df_regions["chr"] = (
        df_regions["chr"].astype(str).str.replace("^chr", "", regex=True)
    )

    df_genes["start_int"] = df_genes["start"].astype(int)
    df_genes["end_int"] = df_genes["end"].astype(int)
    df_regions["start_int"] = df_regions["start"].astype(int)
    df_regions["end_int"] = df_regions["end"].astype(int)

    df_regions["seg_size_bp"] = df_regions["end_int"] - df_regions["start_int"] + 1

    df_genes["seg_start"] = np.nan
    df_genes["seg_end"] = np.nan
    df_genes["seg_size_bp"] = np.nan

    # MAF columns (if present)
    maf_cols: list[str] = []
    for col in ["maf.expected", "maf.observed"]:
        if col in df_regions.columns:
            maf_cols.append(col)
            df_genes[col] = np.nan

    # Additional segment-level metrics to carry over (if present)
    # NOTE: 'type' is intentionally NOT included here so we don't get 'LOH' in type.
    seg_extra_cols: list[str] = []
    for col in ["seg.mean", "num.snps", "C", "M", "M.flagged"]:
        if col in df_regions.columns:
            seg_extra_cols.append(col)
            df_genes[col] = np.nan

    for chr_name, genes_chr in df_genes.groupby("chr"):
        regs_chr = df_regions[df_regions["chr"] == chr_name]
        if regs_chr.empty or genes_chr.empty:
            continue

        regs_chr = regs_chr.sort_values("start_int")
        seg_starts = regs_chr["start_int"].to_numpy()
        seg_ends = regs_chr["end_int"].to_numpy()
        seg_sizes = regs_chr["seg_size_bp"].to_numpy()

        maf_arrays: dict[str, np.ndarray] = {}
        for col in maf_cols:
            maf_arrays[col] = regs_chr[col].to_numpy()

        extra_arrays: dict[str, np.ndarray] = {}
        for col in seg_extra_cols:
            extra_arrays[col] = regs_chr[col].to_numpy()

        for idx, row in genes_chr.iterrows():
            g_start = row["start_int"]
            g_end = row["end_int"]

            # segments that overlap the gene
            overlap_mask = (seg_starts <= g_end) & (seg_ends >= g_start)
            if not overlap_mask.any():
                continue

            # compute overlap size with each overlapping segment
            o_starts = np.maximum(seg_starts[overlap_mask], g_start)
            o_ends = np.minimum(seg_ends[overlap_mask], g_end)
            overlaps = o_ends - o_starts + 1

            best_idx_local = overlaps.argmax()
            seg_idx = np.where(overlap_mask)[0][best_idx_local]

            seg_start = seg_starts[seg_idx]
            seg_end = seg_ends[seg_idx]
            seg_size = seg_sizes[seg_idx]

            df_genes.at[idx, "seg_start"] = int(seg_start)
            df_genes.at[idx, "seg_end"] = int(seg_end)
            df_genes.at[idx, "seg_size_bp"] = int(seg_size)

            for col in maf_cols:
                df_genes.at[idx, col] = maf_arrays[col][seg_idx]

            for col in seg_extra_cols:
                df_genes.at[idx, col] = extra_arrays[col][seg_idx]

    return df_genes


# =============================================================================
# PON spread & weight significance per gene
# =============================================================================


def add_pon_spread_significance(
    df_genes: pd.DataFrame,
    df_cnr: pd.DataFrame,
    df_pon: pd.DataFrame,
    gene_col_genes: str = "gene.symbol",
    gene_col_cnr: str = "gene",
    chr_col: str = "chromosome",
    spread_col: str = "spread",
    noise_factor_borderline: float = 1.0,
    noise_factor_strong: float = 2.0,
) -> pd.DataFrame:
    """
    For each gene in df_genes, compute PON-based noise and CNVkit weight summaries
    from df_cnr + df_pon and add:

      - pon_spread_median
      - mean_vs_spread
      - pon_spread_flag  (within_noise / borderline / beyond_pon_noise / unknown)
    """
    df_genes = df_genes.copy()
    if df_pon.empty:
        for col in ["pon_spread_median", "mean_vs_spread", "pon_spread_flag"]:
            if col not in df_genes.columns:
                df_genes[col] = pd.NA
        return df_genes

    merged_bins = merge_cnr_pon(df_cnr, df_pon, chr_col=chr_col, spread_col=spread_col)

    merged_bins = merged_bins[
        merged_bins[gene_col_cnr].notna()
        & (merged_bins[gene_col_cnr] != "Antitarget")
        & (merged_bins[gene_col_cnr] != "-")
    ].copy()

    merged_bins[gene_col_cnr] = merged_bins[gene_col_cnr].astype(str)

    merged_expanded = merged_bins.assign(
        **{
            gene_col_cnr: merged_bins[gene_col_cnr]
            .str.split(",")
            .apply(lambda lst: [g.strip() for g in lst if g.strip()])
        }
    ).explode(gene_col_cnr)

    agg_dict = {
        spread_col: [
            ("pon_spread_median", "median"),
        ]
    }

    agg = (
        merged_expanded.groupby(gene_col_cnr)
        .agg(
            **{
                name: pd.NamedAgg(column=col, aggfunc=func)
                for col, specs in agg_dict.items()
                for name, func in specs
            }
        )
        .reset_index()
    )

    df_genes = df_genes.merge(
        agg,
        left_on=gene_col_genes,
        right_on=gene_col_cnr,
        how="left",
    ).drop(columns=[gene_col_cnr])

    # ---------------------------
    # compute mean_vs_spread metric
    # ---------------------------
    eps = 1e-6
    if "gene.mean" in df_genes.columns:
        denom = df_genes["pon_spread_median"].replace(0, eps) + eps
        df_genes["mean_vs_spread"] = df_genes["gene.mean"].abs() / denom
    else:
        df_genes["mean_vs_spread"] = np.nan

    # ---------------------------
    # CNV noise classification
    # ---------------------------
    def classify_row(row: pd.Series) -> str:
        s = row["mean_vs_spread"]
        if pd.isna(s):
            return "unknown"
        if s < noise_factor_borderline:
            return "within_noise"
        if s < noise_factor_strong:
            return "borderline"
        return "beyond_pon_noise"

    df_genes["pon_spread_flag"] = df_genes.apply(classify_row, axis=1)

    # ---------------------------
    # round all relevant metrics
    # ---------------------------
    round_cols = [
        "pon_spread_median",
        "mean_vs_spread",
    ]
    round_cols = [c for c in round_cols if c in df_genes.columns]
    df_genes[round_cols] = df_genes[round_cols].round(3)

    return df_genes


# =============================================================================
# refGene / exon coverage
# =============================================================================


def load_refgene_exons(
    refgene_path: str | Path,
    transcript_selection: str = "longest_tx",
) -> Dict[Tuple[str, str], dict]:
    """
    Parse refGene 'flat' file and select ONE primary transcript per gene,
    then build strand-aware (5'->3') exon intervals.

    Returns:
        {
          (chrom_no_chr, gene_symbol): {
              "transcript": transcript_id,
              "strand": "+" or "-",
              "exons": [(start, end), ...]  # ordered 5'->3'
          },
          ...
        }

    transcript_selection:
      - "longest_tx": choose transcript with the largest (txEnd - txStart)
    """
    cols = [
        "gene_symbol",  # 0
        "transcript",  # 1
        "chrom",  # 2
        "strand",  # 3
        "txStart",  # 4
        "txEnd",  # 5
        "cdsStart",  # 6
        "cdsEnd",  # 7
        "exonCount",  # 8
        "exonStarts",  # 9
        "exonEnds",  # 10
    ]
    rg = pd.read_csv(refgene_path, sep="\t", header=None, names=cols)

    rg["chrom"] = rg["chrom"].astype(str).str.replace("^chr", "", regex=True)
    rg["txStart"] = rg["txStart"].astype(int)
    rg["txEnd"] = rg["txEnd"].astype(int)

    exon_map: Dict[Tuple[str, str], dict] = {}

    grouped = rg.groupby(["chrom", "gene_symbol"], as_index=False)

    for (chrom, gene), sub in grouped:
        if transcript_selection == "longest_tx":
            sub = sub.copy()
            sub["tx_len"] = sub["txEnd"] - sub["txStart"]
            best = sub.loc[sub["tx_len"].idxmax()]
        else:
            best = sub.iloc[0]

        transcript = best["transcript"]
        strand = best["strand"]

        starts = [int(x) for x in str(best["exonStarts"]).split(",") if x != ""]
        ends = [int(x) for x in str(best["exonEnds"]).split(",") if x != ""]
        exons: List[Tuple[int, int]] = list(zip(starts, ends))

        if not exons:
            continue

        exons.sort(key=lambda x: x[0])

        if strand == "-":
            exons = exons[::-1]

        exon_map[(chrom, gene)] = {
            "transcript": transcript,
            "strand": strand,
            "exons": exons,
        }

    return exon_map


def annotate_genes_with_exons(
    df_genes: pd.DataFrame,
    exon_map: Dict[Tuple[str, str], dict],
    coverage_threshold: float = 0.95,
) -> pd.DataFrame:
    """
    For each gene in df_genes, use the segment (seg_start, seg_end) to annotate:
      - transcript: the transcript used for exon numbering
      - exon_coverage: fraction of total exonic bases covered by the CNV segment
      - exons_hit: description string ("whole gene",
                    "exons 3–5 of 10", "exons 2,4,7 of 10", etc.)

    Exon numbering is strand-aware (biological 5'->3').
    """
    df_genes = df_genes.copy()

    if "chr" not in df_genes.columns:
        raise ValueError("df_genes must have a 'chr' column.")
    if "gene.symbol" not in df_genes.columns and "gene" not in df_genes.columns:
        raise ValueError("df_genes must have 'gene.symbol' or 'gene' column.")
    if "seg_start" not in df_genes.columns or "seg_end" not in df_genes.columns:
        raise ValueError("df_genes must have 'seg_start' and 'seg_end' columns.")

    df_genes["chr"] = df_genes["chr"].astype(str).str.replace("^chr", "", regex=True)

    df_genes["transcript"] = ""
    df_genes["exon_coverage"] = np.nan
    df_genes["exons_hit"] = ""

    gene_col = "gene.symbol" if "gene.symbol" in df_genes.columns else "gene"

    for idx, row in df_genes.iterrows():
        chrom = str(row["chr"])
        gene = row[gene_col]

        entry = exon_map.get((chrom, gene))
        if entry is None:
            continue

        exons = entry["exons"]
        transcript = entry["transcript"]

        if pd.isna(row["seg_start"]) or pd.isna(row["seg_end"]):
            continue

        seg_start = int(row["seg_start"])
        seg_end = int(row["seg_end"])

        total_exonic_bp = 0
        covered_exonic_bp = 0
        hit_exon_indices = []

        for i, (e_start, e_end) in enumerate(exons, start=1):
            exon_len = e_end - e_start + 1
            total_exonic_bp += exon_len

            ov_start = max(e_start, seg_start)
            ov_end = min(e_end, seg_end)
            if ov_end >= ov_start:
                covered_exonic_bp += ov_end - ov_start + 1
                hit_exon_indices.append(i)

        if total_exonic_bp == 0:
            continue

        cov = covered_exonic_bp / total_exonic_bp

        df_genes.at[idx, "transcript"] = transcript
        df_genes.at[idx, "exon_coverage"] = round(cov, 3)

        if cov >= coverage_threshold:
            desc = "whole gene"
        elif not hit_exon_indices:
            desc = "no coding exons"
        else:
            N = len(exons)
            if len(hit_exon_indices) > 1 and max(hit_exon_indices) - min(
                hit_exon_indices
            ) + 1 == len(hit_exon_indices):
                desc = f"exons {min(hit_exon_indices)}–{max(hit_exon_indices)} of {N}"
            else:
                desc = (
                    "exons " + ",".join(str(i) for i in hit_exon_indices) + f" of {N}"
                )

        df_genes.at[idx, "exons_hit"] = desc

    return df_genes


# =============================================================================
# Final gene table builder (CNR base + optional PureCN + PON)
# =============================================================================


def build_gene_table(
    loh_regions_path: str | Path | None,
    loh_genes_path: str | Path | None,
    cnr_path: str | Path,
    pon_path: str | Path | None,
    refgene_path: str | Path,
    cytoband_path: str | Path,
) -> pd.DataFrame:
    """
    Construction of the per-gene CNV table for reporting:

      - Base: CNR-derived gene stats (always available)
      - Optional: LOH segment overlap from PureCN regions (seg_start, seg_end, seg_size_bp, maf.*, seg.mean, C/M/etc.)
      - Optional: gene-level LOH/type/C/M from PureCN genes file
      - Optional: PON spread & weights (pon_spread_median, mean_vs_spread, etc.)
      - exon coverage and transcript
      - cytoband
      - cleaned 'type' column (AMPLIFICATION/DELETION where missing)
      - rounded / cleaned numeric columns
    """
    # 1) base gene table from CNR
    df_genes = build_genes_from_cnr(cnr_path)
    if df_genes.empty:
        return df_genes

    # ensure seg_* columns exist even if we never annotate segments
    for col in ["seg_start", "seg_end", "seg_size_bp"]:
        if col not in df_genes.columns:
            df_genes[col] = np.nan

    df_cnr = safe_read_csv(cnr_path, sep="\t")
    df_pon = (
        safe_read_csv(pon_path, sep="\t")
        if pon_path is not None and Path(pon_path).is_file()
        else pd.DataFrame()
    )

    # 2) overlay PureCN LOH segments (if available)
    if loh_regions_path is not None and Path(loh_regions_path).is_file():
        df_regions = safe_read_csv(loh_regions_path, sep=",")
        if not df_regions.empty:
            df_genes = annotate_genes_with_segments(df_genes, df_regions)

    # drop helper columns we don't want in final table (if they sneaked in)
    drop_cols = [
        "Sampleid",
        "C.flagged",
        "seg.id",
        "start_int",
        "end_int",
        "breakpoints",
    ]
    drop_cols = [c for c in drop_cols if c in df_genes.columns]
    df_genes = df_genes.drop(columns=drop_cols)

    # 2b) overlay PureCN gene-level LOH (and optionally other columns)
    if loh_genes_path is not None and Path(loh_genes_path).is_file():
        df_loh_genes = safe_read_csv(loh_genes_path, sep=",")
        if not df_loh_genes.empty:
            # normalize chr
            if "chr" in df_loh_genes.columns:
                df_loh_genes["chr"] = (
                    df_loh_genes["chr"].astype(str).str.replace("^chr", "", regex=True)
                )
            elif "chromosome" in df_loh_genes.columns:
                df_loh_genes["chr"] = (
                    df_loh_genes["chromosome"]
                    .astype(str)
                    .str.replace("^chr", "", regex=True)
                )

            # make sure gene.symbol exists
            gene_col = (
                "gene.symbol" if "gene.symbol" in df_loh_genes.columns else "gene"
            )
            df_loh_genes = df_loh_genes.rename(columns={gene_col: "gene.symbol"})

            # LOH per (chr, gene.symbol): any(TRUE)
            if "loh" in df_loh_genes.columns:
                tmp = df_loh_genes.copy()
                tmp["loh_str"] = tmp["loh"].astype(str).str.upper()
                tmp["loh_flag"] = tmp["loh_str"] == "TRUE"

                gene_loh = tmp.groupby(["chr", "gene.symbol"], as_index=False)[
                    "loh_flag"
                ].any()
                gene_loh["loh"] = gene_loh["loh_flag"].replace(
                    {True: True, False: False}
                )
                gene_loh = gene_loh[["chr", "gene.symbol", "loh"]]

                # Merge into df_genes; this becomes canonical LOH
                df_genes = df_genes.merge(
                    gene_loh,
                    on=["chr", "gene.symbol"],
                    how="left",
                )

            # optional: other gene-level columns (first per group)
            extra_cols = []
            for col in ["M.flagged", "type", "C", "M"]:
                if col in df_loh_genes.columns:
                    extra_cols.append(col)

            if extra_cols:
                gene_extra = df_loh_genes.groupby(
                    ["chr", "gene.symbol"], as_index=False
                )[extra_cols].first()
                df_genes = df_genes.merge(
                    gene_extra,
                    on=["chr", "gene.symbol"],
                    how="left",
                    suffixes=("", "_purecn"),
                )
    else:
        if "loh" not in df_genes.columns:
            df_genes["loh"] = pd.NA

    # 3) PON spread + weights
    df_genes = add_pon_spread_significance(df_genes, df_cnr, df_pon)

    # 4) exon coverage
    exon_map = load_refgene_exons(refgene_path)
    df_genes = annotate_genes_with_exons(df_genes, exon_map, coverage_threshold=0.95)

    # 5) type cleanup (keep only AMP/DEL-like; LOH stays in its own column)
    if "type" in df_genes.columns:
        df_genes["type"] = df_genes["type"].astype("string")

        # treat obvious missing markers as NA
        df_genes["type"] = df_genes["type"].replace(
            ["NA", "NaN", "nan", "None", ".", ""],
            pd.NA,
        )

        # any LOH-like type string should NOT live in 'type' – we have a separate loh column
        type_str = df_genes["type"].astype(str).str.upper()
        loh_like = type_str.str.contains("LOH", na=False)
        df_genes.loc[loh_like, "type"] = pd.NA

        # also, if loh column says TRUE, force type to NA (LOH is handled separately)
        if "loh" in df_genes.columns:
            loh_true = df_genes["loh"].astype(str).str.upper() == "TRUE"
            df_genes.loc[loh_true, "type"] = pd.NA

        # now fill missing type from C where possible
        if "C" in df_genes.columns:
            missing_type = df_genes["type"].isna()
            gain_mask = missing_type & df_genes["C"].notna() & (df_genes["C"] > 2)
            loss_mask = missing_type & df_genes["C"].notna() & (df_genes["C"] < 2)

            df_genes.loc[gain_mask, "type"] = "AMPLIFICATION"
            df_genes.loc[loss_mask, "type"] = "DELETION"

    # 6) cytoband
    cyto = load_cytobands(cytoband_path)
    df_genes = annotate_genes_with_cytoband(df_genes, cyto)

    # 6b) sort genes by chromosome and start
    if "chr" in df_genes.columns and "start" in df_genes.columns:

        def _chr_order(v: str) -> int:
            s = str(v).replace("chr", "")
            s = s.upper()
            if s.isdigit():
                return int(s)
            if s == "X":
                return 23
            if s == "Y":
                return 24
            return 25  # anything weird goes last

        df_genes["chr_sort"] = df_genes["chr"].astype(str).map(_chr_order)
        # use seg_start if you prefer segment start; fall back to start
        sort_pos_col = "seg_start" if "seg_start" in df_genes.columns else "start"
        # make sure we sort numerically even if seg_start is currently string/NA
        df_genes["_sort_pos"] = pd.to_numeric(df_genes[sort_pos_col], errors="coerce")
        df_genes = df_genes.sort_values(["chr_sort", "_sort_pos"], kind="mergesort")
        df_genes = df_genes.drop(columns=["chr_sort", "_sort_pos"])

    # 7) column order
    cols = [
        "gene.symbol",
        "transcript",
        "number.targets",
        "chr",
        "start",
        "end",
        "seg_start",
        "seg_end",
        "seg_size_bp",
        "cytoband",
        "exon_coverage",
        "exons_hit",
        "type",
        "loh",
        "C",
        "M",
        "M.flagged",
        "seg.mean",
        "num.snps",
        "maf.expected",
        "maf.observed",
        "gene.mean",
        "gene.min",
        "gene.max",
        "pon_spread_median",
        "weight.mean",
        "mean_vs_spread",
        "pon_spread_flag",
    ]
    cols = [c for c in cols if c in df_genes.columns]
    df_genes = df_genes[cols]

    # remove pseudo-genes (if any slipped in)
    if "gene.symbol" in df_genes.columns:
        df_genes = df_genes[~df_genes["gene.symbol"].isin(["Antitarget", "-"])]

    # rounding for some numeric columns (3 decimals)
    round_cols = [
        "maf.expected",
        "maf.observed",
        "gene.mean",
        "gene.min",
        "gene.max",
        "seg.mean",
    ]
    round_cols = [c for c in round_cols if c in df_genes.columns]
    if round_cols:
        df_genes.loc[:, round_cols] = df_genes.loc[:, round_cols].round(3)

    # pretty integers for seg_start/seg_end/size
    int_cols = ["seg_start", "seg_end", "seg_size_bp"]
    int_cols = [c for c in int_cols if c in df_genes.columns]
    if int_cols:
        df_genes[int_cols] = (
            df_genes[int_cols]
            .apply(pd.to_numeric, errors="coerce")
            .astype("Int64")
            .astype(str)
        )

    return df_genes
