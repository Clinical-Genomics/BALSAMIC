import math
from pathlib import Path
from typing import Dict, Tuple, List, Optional, Set
from pandas.errors import EmptyDataError
import fitz
import matplotlib.pyplot as plt
from matplotlib import colormaps
import numpy as np
import pandas as pd
from BALSAMIC.constants.analysis import Gender
import re
from collections.abc import Collection

# =============================================================================
# Generic helpers
# =============================================================================


def pdf_first_page_to_png(pdf_path: str, png_path: str, dpi: int = 300) -> None:
    doc = fitz.open(pdf_path)
    page = doc[0]
    page.get_pixmap(dpi=dpi).save(png_path)
    doc.close()

def load_cancer_gene_set(path: str | Path,
                         min_occurrence: int = 1,
                         only_annotated: bool = True) -> set[str]:
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
    agg = (
        cnr.groupby([chr_col, "gene.symbol"], as_index=False)
        .agg(
            start=("start", "min"),
            end=("end", "max"),
            number_targets=("log2", "count"),
            gene_mean=("log2", "mean"),
            gene_min=("log2", "min"),
            gene_max=("log2", "max"),
            weight_mean=("weight", "mean"),
        )
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
    cnr_path: Path,
    vcf_path: Path,
    gene_seg_df: pd.DataFrame,
    outdir: Path,
    case_id: str,
    pon_path: Optional[Path] = None,
    include_y: bool = False,
    weight_thresh: float = 0.1,
    pct_spread: float = 0.95,
    window: int = 5,
    anti_factor: float = 0.15,
    base_label_offset: float = 1.5,
    neutral_target_factor: float = 0.3,
    highlight_only_cancer: bool = False,
    y_abs_max: float = 3.0,
    focus_genes: Optional[Collection[str]] = None,
    focus_padding_bp: int = 0,
) -> None:
    """
    Create one PNG per chromosome with:

      - PON spread band (if PON available)
      - smoothed log2 coverage (CNR)
      - segments reconstructed from `gene_seg_df`
      - CNV / LOH / PON-significant genes coloured & labelled
      - BAF panel under each chromosome (from VCF)

    Additional logic:

      - Only genes with >= 3 targets are highlighted.
      - If a gene is highlighted *only* because it is PON-significant,
        it must also have is_cancer_gene == True.
      - If highlight_only_cancer is True, then all highlighted genes
        must be cancer genes (for CNV/LOH as well).
      - If focus_genes is provided, only bins / segments overlapping the
        genomic span of those genes are shown (per chromosome).
    """

    MIN_GENE_TARGETS = 3
    outdir.mkdir(parents=True, exist_ok=True)

    # Normalise focus_genes
    if focus_genes is not None:
        focus_genes_set = {str(g) for g in focus_genes}
        zoom_tag = "genes_" + "_".join(sorted(focus_genes_set))
    else:
        focus_genes_set = None
        zoom_tag = None

    # ---------------------------------
    # 0. Chromosome order and VCF
    # ---------------------------------
    chr_order = [str(i) for i in range(1, 23)] + ["X"]
    if include_y:
        chr_order.append("Y")

    # normalize gene_seg_df chr column
    gdf = gene_seg_df.copy()
    if "chr" not in gdf.columns:
        raise ValueError("gene_seg_df must contain a 'chr' column.")
    gdf["chr"] = gdf["chr"].astype(str).str.replace("^chr", "", regex=True)
    gdf = gdf[gdf["chr"].isin(chr_order)]

    # Normalise targets column name if present
    if "n_targets" in gdf.columns:
        targets_col = "n_targets"
    elif "n.targets" in gdf.columns:
        gdf = gdf.rename(columns={"n.targets": "n_targets"})
        targets_col = "n_targets"
    else:
        targets_col = None

    # Normalise cancer-gene column to a boolean if present
    if "is_cancer_gene" in gdf.columns:
        gdf["is_cancer_gene_bool"] = gdf["is_cancer_gene"].fillna(False).astype(bool)
    else:
        gdf["is_cancer_gene_bool"] = False

    # Load VCF with BAF
    vcf = load_vcf_with_baf(
        vcf_path=vcf_path,
        chr_order=chr_order,
    )

    # ---------------------------------
    # 1. Load CNR and optional PON
    # ---------------------------------
    cnr = safe_read_csv(cnr_path, sep="\t")
    if cnr.empty:
        return

    chr_col = "chromosome" if "chromosome" in cnr.columns else "chr"
    cnr[chr_col] = cnr[chr_col].astype(str).str.replace("^chr", "", regex=True)
    cnr = cnr[cnr[chr_col].isin(chr_order)]

    if not {"start", "end", "log2"}.issubset(cnr.columns):
        missing = {"start", "end", "log2"} - set(cnr.columns)
        raise ValueError(f"CNR missing required columns: {missing}")

    if "gene" not in cnr.columns:
        cnr["gene"] = ""

    use_pon = False
    merged = None

    if pon_path is not None and Path(pon_path).is_file():
        pon = safe_read_csv(pon_path, sep="\t")
        if not pon.empty and "spread" in pon.columns:
            try:
                merged = merge_cnr_pon(
                    cnr,
                    pon,
                    chr_col=chr_col,
                    spread_col="spread",
                )
                if not merged.empty:
                    use_pon = True
            except Exception:
                use_pon = False

    if not use_pon:
        merged = cnr.copy()
        merged[chr_col] = merged[chr_col].astype(str).str.replace("^chr", "", regex=True)
        if "spread" not in merged.columns:
            merged["spread"] = np.nan

    merged = merged[merged[chr_col].isin(chr_order)]

    # Global PON spread threshold (used to drop noisy bins)
    if use_pon:
        spread_thresh = merged["spread"].quantile(pct_spread)
    else:
        spread_thresh = np.inf

    # ---------------------------------
    # 2. Chunk-level CNV classification in gene_seg_df
    # ---------------------------------
    def _is_loh_or_cnv_chunk(row: pd.Series) -> bool:
        if "loh_flag" in row.index and pd.notna(row["loh_flag"]):
            if str(row["loh_flag"]).strip().upper() == "TRUE":
                return True
        for col in ("cnvkit_cnv_call", "purecn_cnv_call"):
            if col in row.index and pd.notna(row[col]):
                val = str(row[col]).strip().upper()
                if val in ("DELETION", "AMPLIFICATION"):
                    return True
        return False

    def _is_pon_signif_chunk(row: pd.Series) -> bool:
        if "pon_call" in row.index and pd.notna(row["pon_call"]):
            return str(row["pon_call"]).strip().lower() == "significant"
        return False

    gdf["is_loh_or_cnv"] = gdf.apply(_is_loh_or_cnv_chunk, axis=1)
    gdf["is_pon_signif"] = gdf.apply(_is_pon_signif_chunk, axis=1)
    gdf["cnv_flag"] = gdf["is_loh_or_cnv"] | gdf["is_pon_signif"]

    # ---------------------------------
    # 2b. Gene-level aggregation
    # ---------------------------------
    if "gene.symbol" in gdf.columns:
        group_cols = ["chr", "gene.symbol"]
        grouped = gdf.groupby(group_cols, as_index=False)

        agg_dict = {
            "has_loh_or_cnv": ("is_loh_or_cnv", "any"),
            "has_pon_sig": ("is_pon_signif", "any"),
            "is_cancer_gene": ("is_cancer_gene_bool", "any"),
        }
        if targets_col is not None:
            agg_dict["total_targets"] = (targets_col, "sum")

        gene_level = grouped.agg(**agg_dict)

        if "total_targets" not in gene_level.columns:
            gene_level["total_targets"] = np.nan
        gene_level["total_targets"] = gene_level["total_targets"].fillna(0)
        gene_level["is_cancer_gene"] = (
            gene_level["is_cancer_gene"].fillna(False).astype(bool)
        )

        mask_base = gene_level["has_loh_or_cnv"] | gene_level["has_pon_sig"]
        mask_base &= gene_level["total_targets"] >= MIN_GENE_TARGETS

        pon_only = gene_level["has_pon_sig"] & ~gene_level["has_loh_or_cnv"]
        mask_pon_cancer = ~pon_only | (pon_only & gene_level["is_cancer_gene"])

        mask_highlight = mask_base & mask_pon_cancer

        if highlight_only_cancer:
            mask_highlight &= gene_level["is_cancer_gene"]

        gene_level["highlight_gene"] = mask_highlight
    else:
        gene_level = None

    # ---------------------------------
    # 3. Per-chromosome plotting
    # ---------------------------------
    y_clip = float(y_abs_max)
    y_lim_chr = (-y_clip, y_clip)

    for chr_name in chr_order:
        sub = merged[merged[chr_col] == chr_name].copy()
        if sub.empty:
            continue

        sub = sub.sort_values("start")

        sub["type"] = np.where(sub["gene"] == "Antitarget", "Antitarget", "Target")

        if "weight" in sub.columns:
            w_mask = sub["weight"] > weight_thresh
        else:
            w_mask = np.ones(len(sub), dtype=bool)

        if use_pon:
            s_mask = sub["spread"] <= spread_thresh
        else:
            s_mask = np.ones(len(sub), dtype=bool)

        sub = sub[w_mask & s_mask]
        if sub.empty:
            continue

        g_chr = gdf[gdf["chr"] == chr_name].copy()
        if g_chr.empty:
            g_chr = pd.DataFrame(columns=gdf.columns)

        # ---------- NEW: restrict to focus_genes window, if requested ----------
        if focus_genes_set is not None:
            g_focus_chr = g_chr[g_chr["gene.symbol"].isin(focus_genes_set)].copy()
            if g_focus_chr.empty:
                # No genes of interest on this chromosome → skip
                continue

            # Try to use region_start/region_end if available, else seg_start/seg_end
            if {"region_start", "region_end"}.issubset(g_focus_chr.columns):
                start_min = g_focus_chr["region_start"].min()
                end_max = g_focus_chr["region_end"].max()
            elif {"seg_start", "seg_end"}.issubset(g_focus_chr.columns):
                start_min = g_focus_chr["seg_start"].min()
                end_max = g_focus_chr["seg_end"].max()
            else:
                # Fallback: use bin-level positions from CNR
                bins_focus = sub[sub["gene"].isin(focus_genes_set)]
                if bins_focus.empty:
                    continue
                start_min = bins_focus["start"].min()
                end_max = bins_focus["end"].max()

            start_min = max(0, int(start_min) - int(focus_padding_bp))
            end_max = int(end_max) + int(focus_padding_bp)

            # Restrict CNR bins to this genomic window
            sub = sub[(sub["end"] >= start_min) & (sub["start"] <= end_max)]
            if sub.empty:
                continue

            # Restrict gene chunks to genes of interest (this chromosome)
            g_chr = g_focus_chr

            # BAF will be filtered later after we compute x-coords, but we
            # already know the genomic span for filtering POS.
        # ----------------------------------------------------------------------

        # Determine which genes to highlight on this chromosome
        if gene_level is not None:
            if focus_genes_set is not None:
                gene_level_chr = gene_level[
                    (gene_level["chr"] == chr_name)
                    & (gene_level["gene.symbol"].isin(focus_genes_set))
                ]
            else:
                gene_level_chr = gene_level[gene_level["chr"] == chr_name]

            g_chr_gene = gene_level_chr[gene_level_chr["highlight_gene"]].copy()
            highlighted_genes = g_chr_gene["gene.symbol"].astype(str).unique()
        else:
            if highlight_only_cancer and "is_cancer_gene" in g_chr.columns:
                mask_gene = g_chr["cnv_flag"] & g_chr["is_cancer_gene"].fillna(False)
            else:
                mask_gene = g_chr["cnv_flag"]
            g_chr_cnv = g_chr[mask_gene].copy()
            highlighted_genes = (
                g_chr_cnv["gene.symbol"].dropna().astype(str).unique()
                if "gene.symbol" in g_chr_cnv.columns
                else np.array([])
            )

        if focus_genes_set is not None:
            # If we zoom to genes, always treat those genes as "visible", even
            # if not flagged as CNV – you can adjust this if you want.
            highlighted_genes = np.union1d(
                highlighted_genes, list(focus_genes_set)
            )

        if highlighted_genes.size > 0 and "gene.symbol" in g_chr.columns:
            g_chr_high = g_chr[g_chr["gene.symbol"].isin(highlighted_genes)].copy()
        else:
            g_chr_high = pd.DataFrame(columns=g_chr.columns)

        # ---------- Segment info ----------
        segs_chr = pd.DataFrame()
        if {"seg_start", "seg_end"}.issubset(g_chr.columns):
            segs_chr = (
                g_chr.dropna(subset=["seg_start", "seg_end"])
                .groupby(["chr", "seg_start", "seg_end"], as_index=False)
                .agg(
                    seg_log2=("seg_log2", "first")
                    if "seg_log2" in g_chr.columns
                    else ("region_start", "first"),
                    seg_C=("seg_cn", "first")
                    if "seg_cn" in g_chr.columns
                    else ("region_start", "first"),
                    loh_seg_mean=("loh_seg_mean", "first")
                    if "loh_seg_mean" in g_chr.columns
                    else ("region_start", "first"),
                    loh_C=("loh_C", "first")
                    if "loh_C" in g_chr.columns
                    else ("region_start", "first"),
                )
            )
            segs_chr = segs_chr[segs_chr["chr"] == chr_name].sort_values("seg_start")

        # ---------- Variable-width x ----------
        def _bin_width(row: pd.Series) -> float:
            if row["type"] == "Antitarget":
                return anti_factor
            if row["gene"] in highlighted_genes:
                return 1.0
            return neutral_target_factor

        sub["bin_width"] = sub.apply(_bin_width, axis=1)
        sub["x_coord"] = sub["bin_width"].cumsum() - sub["bin_width"] / 2

        # smoothing
        sub = sub.sort_values("x_coord")
        sub["log2_smooth"] = sub["log2"].rolling(window=window, center=True).median()
        if use_pon:
            sub["spread_smooth"] = sub["spread"].rolling(window=window, center=True).median()
            if "pon_log2" in sub.columns:
                sub["pon_log2_smooth"] = sub["pon_log2"].rolling(
                    window=window, center=True
                ).median()
            else:
                sub["pon_log2_smooth"] = 0.0
        else:
            sub["spread_smooth"] = np.nan
            sub["pon_log2_smooth"] = 0.0

        sub["log2_clipped"] = sub["log2"].clip(-y_clip, y_clip)
        sub["log2_smooth_clipped"] = sub["log2_smooth"].clip(-y_clip, y_clip)

        bin_starts = sub["start"].values
        x_coords = sub["x_coord"].values

        def pos_to_xcoord(pos: int) -> float:
            idx = np.searchsorted(bin_starts, pos, side="right") - 1
            if idx < 0:
                idx = 0
            if idx >= len(x_coords):
                idx = len(x_coords) - 1
            return x_coords[idx]

        # BAF for this chromosome, optionally restricted to focus window
        baf_chr = vcf[vcf["CHROM"] == chr_name].copy()
        baf_chr = baf_chr.sort_values("POS")
        if focus_genes_set is not None and not baf_chr.empty:
            # We know our current CNR span from sub
            start_span = sub["start"].min()
            end_span = sub["end"].max()
            baf_chr = baf_chr[
                (baf_chr["POS"] >= start_span) & (baf_chr["POS"] <= end_span)
            ]
        if not baf_chr.empty:
            baf_chr["x_coord"] = baf_chr["POS"].apply(pos_to_xcoord)

        # gene positions
        gene_to_color: dict[str, tuple] = {}
        genes_chr = pd.DataFrame()

        if highlighted_genes.size > 0 and not g_chr_high.empty:
            if {"gene.symbol", "region_start"}.issubset(g_chr_high.columns):
                gene_pos = (
                    g_chr_high.groupby("gene.symbol", as_index=False)["region_start"]
                    .min()
                    .rename(columns={"region_start": "gene_start"})
                )
            else:
                gene_pos = (
                    g_chr_high.groupby("gene.symbol", as_index=False)["seg_start"]
                    .min()
                    .rename(columns={"seg_start": "gene_start"})
                )

            gene_pos["x_coord"] = gene_pos["gene_start"].apply(pos_to_xcoord)
            genes_chr = gene_pos

            n_genes = len(highlighted_genes)
            cmap = colormaps.get_cmap("tab20").resampled(max(n_genes, 1))
            for i, gname in enumerate(highlighted_genes):
                gene_to_color[gname] = cmap(i)

        # segments x-coordinates
        if not segs_chr.empty:
            segs_chr = segs_chr.copy()
            segs_chr["x_start"] = segs_chr["seg_start"].apply(pos_to_xcoord)
            segs_chr["x_end"] = segs_chr["seg_end"].apply(pos_to_xcoord)

            if "seg_log2" in segs_chr.columns:
                segs_chr["seg_log2_clipped"] = segs_chr["seg_log2"].clip(-y_clip, y_clip)
            if "loh_seg_mean" in segs_chr.columns:
                segs_chr["loh_seg_mean_clipped"] = segs_chr["loh_seg_mean"].clip(
                    -y_clip, y_clip
                )

        # ---------- figure ----------
        fig, (ax1, ax2) = plt.subplots(
            2,
            1,
            figsize=(14, 6),
            sharex=True,
            gridspec_kw={"height_ratios": [2, 1]},
        )

        x = sub["x_coord"]

        if highlighted_genes.size > 0:
            non_cnv = sub[~sub["gene"].isin(highlighted_genes)]
        else:
            non_cnv = sub

        bg_targets = non_cnv[non_cnv["type"] == "Target"]
        bg_antis = non_cnv[non_cnv["type"] == "Antitarget"]

        if not bg_antis.empty:
            ax1.scatter(
                bg_antis["x_coord"],
                bg_antis["log2_clipped"],
                s=3,
                alpha=0.38,
                color="lightgrey",
                label="Antitarget bins",
            )

        if not bg_targets.empty:
            ax1.scatter(
                bg_targets["x_coord"],
                bg_targets["log2_clipped"],
                s=4,
                alpha=0.5,
                color="tab:blue",
                label="Target bins (no highlighted CNV gene)",
            )

        if highlighted_genes.size > 0:
            for gene in highlighted_genes:
                gsub = sub[sub["gene"] == gene]
                if gsub.empty:
                    continue
                color = gene_to_color.get(gene, "black")
                ax1.scatter(
                    gsub["x_coord"],
                    gsub["log2_clipped"],
                    s=8,
                    alpha=0.9,
                    color=color,
                )

        # PON band centred at PON mean
        if (
            use_pon
            and "spread_smooth" in sub.columns
            and "pon_log2_smooth" in sub.columns
            and sub["spread_smooth"].notna().any()
        ):
            band_center = sub["pon_log2_smooth"].fillna(0.0)
            band_half = sub["spread_smooth"].fillna(0.0)
            band_bottom = (band_center - band_half).clip(-y_clip, y_clip)
            band_top = (band_center + band_half).clip(-y_clip, y_clip)

            ax1.fill_between(
                x,
                band_bottom,
                band_top,
                alpha=0.4,
                step="mid",
                label="PON noise band",
            )

        ax1.plot(
            x,
            sub["log2_smooth_clipped"],
            linewidth=1.5,
            alpha=0.9,
            color="tab:green",
            label=f"log2 (median {window} bins)",
        )

        if not segs_chr.empty:
            for _, srow in segs_chr.iterrows():
                xs = srow["x_start"]
                xe = srow["x_end"]

                if (
                    "loh_seg_mean_clipped" in srow.index
                    and pd.notna(srow["loh_seg_mean_clipped"])
                ):
                    y_seg = srow["loh_seg_mean_clipped"]
                elif (
                    "seg_log2_clipped" in srow.index
                    and pd.notna(srow["seg_log2_clipped"])
                ):
                    y_seg = srow["seg_log2_clipped"]
                else:
                    continue

                C_val = np.nan
                if "loh_C" in srow.index and pd.notna(srow["loh_C"]):
                    C_val = srow["loh_C"]
                elif "seg_C" in srow.index and pd.notna(srow["seg_C"]):
                    C_val = srow["seg_C"]

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
        ax1.set_ylabel("log2 / PON band")
        title_suffix = "log2 vs PON spread" if use_pon else "log2 (no PON available)"

        title = f"Chr {chr_name} – {title_suffix}: {case_id}"
        if focus_genes_set is not None:
            title += f"  (genes: {', '.join(sorted(focus_genes_set))})"
        ax1.set_title(title + "\n")
        ax1.legend(loc="upper right", fontsize=8)

        # gene labels
        if not genes_chr.empty:
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
                y_val = sub["log2_smooth_clipped"].iloc[idx_near]
                if pd.isna(y_val):
                    y_val = 0.0

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
        ax2.set_xlabel(
            "Pseudo-position (highlighted genes expanded, other bins compressed)"
        )
        ax2.set_ylabel("BAF")

        plt.tight_layout()

        # file name, with optional gene zoom tag
        if zoom_tag is not None:
            out_png = outdir / f"cnv_chr{chr_name}_{zoom_tag}_segments.png"
        else:
            out_png = outdir / f"cnv_chr{chr_name}_segments.png"

        plt.savefig(out_png, dpi=150)
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


def classify_cnv_from_total_cn_sex_aware(
    cn: float | int | None,
    chrom: str | int,
    sex: Gender,
) -> str:
    """
    Sex-aware classifier based on absolute total copy number.

    - Autosomes (1–22): baseline = 2
    - Chr X: baseline = 2 (female), 1 (male)
    - Chr Y: baseline = 1 (male), 0 (female)

    Returns: 'DELETION', 'AMPLIFICATION', or 'NEUTRAL'
    """
    if pd.isna(cn):
        return ""

    # Clean chromosome string
    chrom_str = str(chrom)
    chrom_str = re.sub(r"^chr", "", chrom_str, flags=re.IGNORECASE)

    # Baseline expectation by sex & chromosome
    if chrom_str.isdigit():
        expected = 2  # autosomes
    elif chrom_str.upper() == "X":
        expected = 2 if sex == Gender.FEMALE else 1
    elif chrom_str.upper() == "Y":
        expected = 1 if sex == Gender.MALE else 0
    else:
        # unknown contig: fall back to 2
        expected = 2

    cn_val = float(cn)
    # Round both to integers for robustness
    cn_int = int(round(cn_val))
    exp_int = int(round(expected))

    if cn_int < exp_int:
        return "DELETION"
    elif cn_int > exp_int:
        return "AMPLIFICATION"
    else:
        return "NEUTRAL"

# =============================================================================
# Final gene table builder (CNR base + optional PureCN + PON)
# =============================================================================

def build_gene_segment_table(
        cnr_path: str | Path,
        cns_path: str | Path,
        cancer_genes: Optional[Set[str]] = None,
        refgene_path: str | Path | None = None,
        transcript_selection: str = "longest_tx",
        pon_path: str | Path | None = None,
        loh_path: str | Path | None = None,
        cytoband_path: str | Path | None = None,
        sex: Gender = None,
) -> pd.DataFrame:
    """
    Build a gene x segment table from CNVkit CNR + CNS (+ optional refGene + optional PON).

    One row ≈ one (gene, segment) chunk:

      - A gene that lies entirely in one CNS segment → 1 row.
      - A gene that spans two CNS segments → 2 rows.
      - Gene parts without any CNS segment (neutral / no-call) form their own row.

    Optional:
      - refgene_path: annotate exons_hit (whole_gene vs exon indices).
      - pon_path: use CNVkit PON .cnn to compute a PON-based deviation per chunk.

    Key columns:

      chr, gene.symbol,
      region_start, region_end,
      n_targets, mean_log2, min_log2, max_log2,
      seg_start, seg_end, seg_log2, seg_baf, seg_cn, seg_cn1, seg_cn2,
      is_cancer_gene,
      exons_hit (if refgene_path),
      pon_mean_log2, pon_mean_spread, pon_effect, pon_z, pon_direction, pon_call (if pon_path).
    """

    # ------------------------- 1. Load CNR ------------------------- #
    cnr = safe_read_csv(cnr_path, sep="\t").copy()
    if cnr.empty:
        return pd.DataFrame()

    cnr_chr_col = "chromosome" if "chromosome" in cnr.columns else "chr"
    cnr[cnr_chr_col] = cnr[cnr_chr_col].astype(str).str.replace("^chr", "", regex=True)

    if not {"start", "end", "log2", "gene"}.issubset(cnr.columns):
        missing = {"start", "end", "log2", "gene"} - set(cnr.columns)
        raise ValueError(f"CNR missing required columns: {missing}")

    has_depth = "depth" in cnr.columns
    has_weight = "weight" in cnr.columns

    # drop pseudo-genes
    cnr = cnr[cnr["gene"].notna()]
    cnr = cnr[~cnr["gene"].isin(["Antitarget", "-"])]

    # explode multi-gene bins
    cnr["gene.symbol"] = (
        cnr["gene"]
        .astype(str)
        .str.split(",")
        .apply(lambda lst: [g.strip() for g in lst if g.strip()])
    )
    cnr = cnr.explode("gene.symbol")

    # normalize column names
    cnr = cnr.rename(columns={cnr_chr_col: "chr"})

    # ------------------------- 1b. Optional: Load PON ------------------------- #
    # Join by (chr, start, end). Assumes PON is built on same bin design.
    if pon_path is not None and Path(pon_path).is_file():
        pon = safe_read_csv(pon_path, sep="\t").copy()
        if not pon.empty:
            pon_chr_col = "chromosome" if "chromosome" in pon.columns else "chr"
            pon[pon_chr_col] = pon[pon_chr_col].astype(str).str.replace("^chr", "", regex=True)
            if not {"start", "end", "log2", "spread"}.issubset(pon.columns):
                raise ValueError("PON file must contain 'start','end','log2','spread' columns.")
            pon = pon.rename(
                columns={
                    pon_chr_col: "chr",
                    "log2": "pon_log2",
                    "spread": "pon_spread",
                }
            )
            cnr = cnr.merge(
                pon[["chr", "start", "end", "pon_log2", "pon_spread"]],
                how="left",
                on=["chr", "start", "end"],
            )

    # ------------------------- 2. Load CNS ------------------------- #
    cns = safe_read_csv(cns_path, sep="\t").copy()
    if cns.empty:
        grouped = (
            cnr.groupby(["chr", "gene.symbol"], as_index=False)
            .agg(
                region_start=("start", "min"),
                region_end=("end", "max"),
                n_targets=("log2", "count"),
                mean_log2=("log2", "mean"),
                min_log2=("log2", "min"),
                max_log2=("log2", "max"),
                mean_depth=("depth", "mean") if has_depth else ("log2", "size"),
                mean_weight=("weight", "mean") if has_weight else ("log2", "size"),
            )
        )
        grouped["segment_id"] = "no_segment"
        grouped["seg_start"] = np.nan
        grouped["seg_end"] = np.nan
        grouped["seg_log2"] = np.nan
        grouped["seg_baf"] = np.nan
        grouped["seg_cn"] = np.nan
        grouped["seg_cn1"] = np.nan
        grouped["seg_cn2"] = np.nan
        grouped["seg_depth"] = np.nan

        if cancer_genes is not None:
            grouped["is_cancer_gene"] = grouped["gene.symbol"].isin(cancer_genes)

        # Simple case: no segments; we won't do exon or PON logic here for now.
        return grouped

    cns_chr_col = "chromosome" if "chromosome" in cns.columns else "chr"
    cns[cns_chr_col] = cns[cns_chr_col].astype(str).str.replace("^chr", "", regex=True)

    if "start" not in cns.columns or "end" not in cns.columns:
        raise ValueError("CNS file must contain 'start' and 'end' columns.")

    base_seg_cols = ["start", "end", "log2", "baf", "cn", "cn1", "cn2", "depth"]
    seg_cols_available = [c for c in base_seg_cols if c in cns.columns]

    cns = cns[[cns_chr_col] + seg_cols_available].copy()
    cns = cns.rename(columns={cns_chr_col: "chr"})
    cns = cns.sort_values(["chr", "start"]).reset_index(drop=True)
    cns["segment_id"] = cns.index.astype(str)

    cns_per_chr: dict[str, pd.DataFrame] = {
        chrom: df_chr.reset_index(drop=True)
        for chrom, df_chr in cns.groupby("chr")
    }

    # -------------------- 3. Annotate bins with segment -------------------- #
    def _assign_segment_to_bins(df_bins: pd.DataFrame, segs: pd.DataFrame) -> pd.DataFrame:
        """
        For each bin, assign the segment that covers the bin *center*.
        If none, seg_* stay NaN and segment_id = "no_segment".
        """
        if segs is None or segs.empty or df_bins.empty:
            df_bins = df_bins.copy()
            df_bins["segment_id"] = "no_segment"
            for col in ["seg_start", "seg_end", "seg_log2", "seg_baf", "seg_cn", "seg_cn1", "seg_cn2", "seg_depth"]:
                df_bins[col] = np.nan
            return df_bins

        seg_starts = segs["start"].values
        seg_ends = segs["end"].values

        seg_id_col = np.full(len(df_bins), "no_segment", dtype=object)
        seg_start_col = np.full(len(df_bins), np.nan)
        seg_end_col = np.full(len(df_bins), np.nan)
        seg_log2_col = np.full(len(df_bins), np.nan)
        seg_baf_col = np.full(len(df_bins), np.nan)
        seg_cn_col = np.full(len(df_bins), np.nan)
        seg_cn1_col = np.full(len(df_bins), np.nan)
        seg_cn2_col = np.full(len(df_bins), np.nan)
        seg_depth_col = np.full(len(df_bins), np.nan)

        centers = ((df_bins["start"].values + df_bins["end"].values) // 2).astype(int)

        for i, center in enumerate(centers):
            mask = (center >= seg_starts) & (center < seg_ends)
            idx = np.flatnonzero(mask)
            if len(idx) == 0:
                continue
            j = idx[0]
            seg = segs.iloc[j]
            seg_id_col[i] = seg["segment_id"]
            seg_start_col[i] = seg["start"]
            seg_end_col[i] = seg["end"]
            if "log2" in segs.columns:
                seg_log2_col[i] = seg["log2"]
            if "baf" in segs.columns:
                seg_baf_col[i] = seg["baf"]
            if "cn" in segs.columns:
                seg_cn_col[i] = seg["cn"]
            if "cn1" in segs.columns:
                seg_cn1_col[i] = seg["cn1"]
            if "cn2" in segs.columns:
                seg_cn2_col[i] = seg["cn2"]
            if "depth" in segs.columns:
                seg_depth_col[i] = seg["depth"]

        df_bins = df_bins.copy()
        df_bins["segment_id"] = seg_id_col
        df_bins["seg_start"] = seg_start_col
        df_bins["seg_end"] = seg_end_col
        df_bins["seg_log2"] = seg_log2_col
        df_bins["seg_baf"] = seg_baf_col
        df_bins["seg_cn"] = seg_cn_col
        df_bins["seg_cn1"] = seg_cn1_col
        df_bins["seg_cn2"] = seg_cn2_col
        df_bins["seg_depth"] = seg_depth_col

        return df_bins

    annotated_bins_list = []
    for chrom, df_chr in cnr.groupby("chr"):
        segs_chr = cns_per_chr.get(chrom)
        annotated = _assign_segment_to_bins(df_chr, segs_chr)
        annotated_bins_list.append(annotated)

    bins = pd.concat(annotated_bins_list, ignore_index=True)

    # -------------------- 4. Collapse to gene x segment -------------------- #
    agg_dict: dict[str, list[str]] = {
        "start": ["min", "max"],  # region_start/region_end
        "log2": ["count", "mean", "min", "max"],
    }
    if has_depth:
        agg_dict["depth"] = ["mean"]
    if has_weight:
        agg_dict["weight"] = ["mean"]
    # Optional PON aggregation
    if "pon_log2" in bins.columns:
        agg_dict["pon_log2"] = ["mean"]
    if "pon_spread" in bins.columns:
        agg_dict["pon_spread"] = ["mean"]

    grouped = (
        bins.groupby(["chr", "gene.symbol", "segment_id"], as_index=False)
        .agg(agg_dict)
    )

    # flatten multi-index columns
    grouped.columns = [
        "_".join(col).strip("_") if isinstance(col, tuple) else col
        for col in grouped.columns
    ]

    grouped = grouped.rename(
        columns={
            "start_min": "region_start",
            "start_max": "region_end",
            "log2_count": "n_targets",
            "log2_mean": "mean_log2",
            "log2_min": "min_log2",
            "log2_max": "max_log2",
            "depth_mean": "mean_depth" if has_depth else "depth_mean",
            "weight_mean": "mean_weight" if has_weight else "weight_mean",
            "pon_log2_mean": "pon_mean_log2",
            "pon_spread_mean": "pon_mean_spread",
        }
    )
    # keep legacy alias if you need it
    if "n_targets" in grouped.columns:
        grouped["n.targets"] = grouped["n_targets"]

    # attach segment attributes (one row per segment_id per chr)
    seg_info = (
        bins.groupby(["chr", "segment_id"], as_index=False)
        .agg(
            seg_start=("seg_start", "first"),
            seg_end=("seg_end", "first"),
            seg_log2=("seg_log2", "first"),
            seg_baf=("seg_baf", "first"),
            seg_cn=("seg_cn", "first"),
            seg_cn1=("seg_cn1", "first"),
            seg_cn2=("seg_cn2", "first"),
            seg_depth=("seg_depth", "first"),
        )
    )

    grouped = grouped.merge(
        seg_info,
        how="left",
        on=["chr", "segment_id"],
    )

    # Cancer gene flag
    if cancer_genes is not None:
        grouped["is_cancer_gene"] = grouped["gene.symbol"].isin(cancer_genes)

    # -------------------- 5. PON-based deviation per chunk -------------------- #
    min_n_for_pon = 2  # ignore PON-driven stat if fewer bins than this

    if "pon_mean_log2" in grouped.columns and "pon_mean_spread" in grouped.columns:
        grouped["pon_effect"] = grouped["mean_log2"] - grouped["pon_mean_log2"]

        def _compute_pon_z(row: pd.Series) -> float:
            eff = row["pon_effect"]
            sigma = row["pon_mean_spread"]
            n = row.get("n_targets", row.get("n.targets", np.nan))

            if pd.isna(eff) or pd.isna(sigma) or sigma <= 0:
                return np.nan
            if pd.isna(n) or n < min_n_for_pon:
                # too few bins for a decent PON-based z
                return np.nan

            # effective sigma of the mean
            sigma_eff = sigma / np.sqrt(float(n))
            if sigma_eff <= 0:
                return np.nan

            return abs(eff) / sigma_eff

        grouped["pon_z"] = grouped.apply(_compute_pon_z, axis=1)

        def _pon_direction(row: pd.Series) -> str:
            eff = row["pon_effect"]
            if pd.isna(eff):
                return ""
            if eff > 0:
                return "gain"
            if eff < 0:
                return "loss"
            return "neutral"

        grouped["pon_direction"] = grouped.apply(_pon_direction, axis=1)

        def _pon_call(row: pd.Series) -> str:
            z = row["pon_z"]
            if pd.isna(z):
                return ""
            if z < 1.5:
                return "noise"
            if z < 3.0:
                return "borderline"
            return "significant"

        grouped["pon_call"] = grouped.apply(_pon_call, axis=1)

    # -------------------- 6. Exon annotation (optional) -------------------- #
    if refgene_path is not None:
        exon_map = load_refgene_exons(
            refgene_path=refgene_path,
            transcript_selection=transcript_selection,
        )

        def _exons_hit(row: pd.Series) -> str:
            key = (str(row["chr"]), str(row["gene.symbol"]))
            info = exon_map.get(key)
            seg_start = row.get("seg_start")
            seg_end = row.get("seg_end")

            if info is None or pd.isna(seg_start) or pd.isna(seg_end):
                return ""

            exons = info["exons"]  # already 5'->3' oriented
            gene_start = min(s for s, e in exons)
            gene_end = max(e for s, e in exons)

            if seg_start <= gene_start and seg_end >= gene_end:
                return "whole_gene"

            hit_indices: list[int] = []
            for idx, (s, e) in enumerate(exons, start=1):
                if e > seg_start and s < seg_end:
                    hit_indices.append(idx)

            if not hit_indices:
                return ""

            # compress contiguous exon numbers: [1,2,3,5] -> "1-3,5"
            ranges = []
            start = prev = hit_indices[0]
            for i in hit_indices[1:]:
                if i == prev + 1:
                    prev = i
                else:
                    ranges.append((start, prev))
                    start = prev = i
            ranges.append((start, prev))

            parts = []
            for a, b in ranges:
                if a == b:
                    parts.append(str(a))
                else:
                    parts.append(f"{a}-{b}")
            return ",".join(parts)

        grouped["exons_hit"] = grouped.apply(_exons_hit, axis=1)

    # -------------------- 7. LOHgenes annotation (optional) -------------------- #
    if loh_path is not None and Path(loh_path).is_file():
        loh = safe_read_csv(loh_path, sep=",").copy()
        if not loh.empty:
            # normalize chromosome naming
            loh_chr_col = "chr" if "chr" in loh.columns else "chromosome"
            loh[loh_chr_col] = loh[loh_chr_col].astype(str).str.replace("^chr", "", regex=True)
            loh = loh.rename(columns={loh_chr_col: "chr"})

            # If multiple samples are present, force user to filter upstream
            if "Sampleid" in loh.columns and loh["Sampleid"].nunique() > 1:
                raise ValueError(
                    "LOHgenes file contains multiple Sampleid values; "
                    "please subset to a single sample before calling build_gene_segment_table."
                )

            # ensure required columns
            required_loh_cols = {"gene.symbol", "chr", "C", "M", "M.flagged", "loh", "seg.mean", "num.snps"}
            missing = required_loh_cols - set(loh.columns)
            if missing:
                raise ValueError(f"LOHgenes file missing columns: {missing}")

            # keep only what we need and rename for clarity
            loh_small = loh[["chr", "gene.symbol", "C", "M", "M.flagged", "loh", "seg.mean", "num.snps"]].copy()
            loh_small = loh_small.rename(
                columns={
                    "C": "loh_C",
                    "M": "loh_M",
                    "M.flagged": "loh_M_flagged",
                    "loh": "loh_flag",
                    "seg.mean": "loh_seg_mean",
                    "num.snps": "loh_num_snps",
                }
            )

            # merge per (chr, gene.symbol); allows multiple segment rows for same gene
            grouped = grouped.merge(
                loh_small,
                how="left",
                on=["chr", "gene.symbol"],
            )

    # 8) cytoband
    cyto = load_cytobands(cytoband_path)
    grouped = annotate_genes_with_cytoband(grouped, cyto)

    # -------------------- 9. Column order & sorting -------------------- #
    cols_order = [
        "chr",
        "region_start",
        "region_end",
        "seg_start",
        "seg_end",
        "cytoband",
        "gene.symbol",
        "n.targets",
        "seg_log2",
        "loh_seg_mean",
        "loh_num_snps",
        "seg_baf",
        "mean_log2",
        "min_log2",
        "max_log2",
        "mean_weight",
        "seg_cn",
        "seg_cn1",
        "seg_cn2",
        "loh_C",
        "loh_M",
        "loh_M_flagged",
        "loh_flag",
        "is_cancer_gene",
        "exons_hit",
        "pon_mean_log2",
        "pon_mean_spread",
        "pon_z",
        "pon_call",
    ]
    cols_order = [c for c in cols_order if c in grouped.columns]
    grouped = grouped[cols_order + (["segment_id"] if "segment_id" in grouped.columns else [])]

    # chromosome sorting: numeric -> X -> Y -> others
    def _chr_key(c):
        try:
            return (0, int(c))
        except ValueError:
            if c in ("X", "x"):
                return (1, 23)
            if c in ("Y", "y"):
                return (1, 24)
            return (2, c)

    grouped["chr_sort"] = grouped["chr"].map(_chr_key)

    grouped = grouped.sort_values(
        by=["chr_sort", "region_start", "region_end"],
        kind="stable",
    ).drop(columns=["chr_sort"])

    # drop helper columns we don't want in final table
    drop_cols = [
        "segment_id",
        "mean_depth",
    ]
    drop_cols = [c for c in drop_cols if c in grouped.columns]
    grouped = grouped.drop(columns=drop_cols)

    # CNVkit-based call (from seg_cn)
    if {"seg_cn", "chr"}.issubset(grouped.columns):
        grouped["cnvkit_cnv_call"] = grouped.apply(
            lambda r: classify_cnv_from_total_cn_sex_aware(
                r["seg_cn"],
                r["chr"],
                sex,
            ),
            axis=1,
        )

    # PureCN-based call (from loh_C)
    if {"loh_C", "chr"}.issubset(grouped.columns):
        grouped["purecn_cnv_call"] = grouped.apply(
            lambda r: classify_cnv_from_total_cn_sex_aware(
                r["loh_C"],
                r["chr"],
                sex,
            ),
            axis=1,
        )

    return grouped
