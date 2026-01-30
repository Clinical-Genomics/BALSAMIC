from __future__ import annotations

# Standard library
import math
from collections.abc import Collection
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

# Third-party
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colormaps
from matplotlib import patheffects as pe

from cnv_report_utils import safe_read_csv


# =============================================================================
# VCF parsing / BAF
# =============================================================================


@dataclass(frozen=True)
class ParsedSample:
    baf: float
    dp: float  # dp as float to allow NaN


def parse_sample_fields(format_str: str, sample_str: str) -> tuple[float, None, float]:
    """
    Parse AD and DP from FORMAT/TUMOR and compute BAF = alt / (ref + alt).

    Returns (BAF, dummy_GT, DP_sample).
    """
    parsed = _parse_sample_fields(format_str, sample_str)
    return parsed.baf, None, parsed.dp


def _parse_sample_fields(format_str: str, sample_str: str) -> ParsedSample:
    """Internal parser returning a structured result."""
    if pd.isna(format_str) or pd.isna(sample_str):
        return ParsedSample(baf=math.nan, dp=math.nan)

    keys = str(format_str).split(":")
    vals = str(sample_str).split(":")
    fmt = dict(zip(keys, vals))

    # DP
    dp_str = fmt.get("DP", None)
    try:
        dp = int(dp_str) if dp_str not in (None, ".", "") else math.nan
    except ValueError:
        dp = math.nan

    # AD
    ad = fmt.get("AD", None)
    if ad is None or ad in (".", ""):
        return ParsedSample(baf=math.nan, dp=dp)

    try:
        counts = [int(x) for x in str(ad).split(",")]
    except ValueError:
        return ParsedSample(baf=math.nan, dp=dp)

    if len(counts) < 2:
        return ParsedSample(baf=math.nan, dp=dp)

    ref_c, alt_c = counts[0], counts[1]  # first alt only
    total = ref_c + alt_c
    if total <= 0:
        return ParsedSample(baf=math.nan, dp=dp)

    return ParsedSample(baf=alt_c / total, dp=dp)


def load_vcf_with_baf(vcf_path: str | Path, chr_order: list[str]) -> pd.DataFrame:
    """
    Load a VCF (single tumor sample column) and compute:
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
    if vcf.empty:
        vcf["BAF"] = []
        vcf["DP_sample"] = []
        return vcf

    vcf["CHROM"] = vcf["CHROM"].astype(str).str.replace("^chr", "", regex=True)
    vcf = vcf[vcf["CHROM"].isin(chr_order)].copy()
    if vcf.empty:
        vcf["BAF"] = []
        vcf["DP_sample"] = []
        return vcf

    parsed = [
        _parse_sample_fields(fmt, sample)
        for fmt, sample in zip(vcf["FORMAT"].tolist(), vcf["TUMOR"].tolist())
    ]
    vcf["BAF"] = [p.baf for p in parsed]
    vcf["DP_sample"] = [p.dp for p in parsed]
    return vcf


# =============================================================================
# CNR + PON merge (dataframe-based)
# =============================================================================


def merge_cnr_pon(
    df_cnr: pd.DataFrame,
    df_pon: pd.DataFrame,
    chr_col: str = "chromosome",
    spread_col: str = "spread",
) -> pd.DataFrame:
    """
    Merge CNVkit tumor CNR with PON on (chr, start, end) and keep only rows
    with non-null PON spread (inner join behavior).
    """
    cnr = df_cnr.copy()
    pon = df_pon.copy()

    if chr_col not in cnr.columns or chr_col not in pon.columns:
        raise ValueError(f"Both cnr and pon must contain chromosome column '{chr_col}'")

    if spread_col not in pon.columns:
        raise ValueError(f"PON is missing required column '{spread_col}'")

    cnr[chr_col] = cnr[chr_col].astype(str).str.replace("^chr", "", regex=True)
    pon[chr_col] = pon[chr_col].astype(str).str.replace("^chr", "", regex=True)

    pon = pon.dropna(subset=[spread_col]).copy()

    merged = pd.merge(
        cnr,
        pon[[chr_col, "start", "end", spread_col]],
        on=[chr_col, "start", "end"],
        how="inner",
        suffixes=("_cnr", "_pon"),
    )

    if spread_col + "_pon" in merged.columns:
        merged[spread_col] = merged[spread_col + "_pon"]

    return merged


# =============================================================================
# Plotting helpers (refactor plot_chromosomes)
# =============================================================================


def _as_chr_order(include_y: bool) -> list[str]:
    """
    Build canonical chromosome ordering list.

    Returns:
      ['1'..'22', 'X'] plus 'Y' if include_y is True.

    Used to standardize filtering and plotting chromosome iteration.
    """
    order = [str(i) for i in range(1, 23)] + ["X"]
    if include_y:
        order.append("Y")
    return order


def _norm_chr_inplace(df: pd.DataFrame, col: str) -> pd.DataFrame:
    """
    Normalize chromosome column to string and strip leading 'chr'.

    Returns a copy of the dataframe with normalized chromosome values.
    Safe to call repeatedly.
    """
    out = df.copy()
    out[col] = out[col].astype(str).str.replace("^chr", "", regex=True)
    return out


def _normalize_gene_tables(
    gene_seg_df: pd.DataFrame,
    chr_order: list[str],
    gene_chunk_df: pd.DataFrame | None,
) -> tuple[pd.DataFrame, pd.DataFrame | None, bool]:
    """
    Normalize gene-level and optional chunk-level gene tables.

    Ensures:
      - 'chr' column exists
      - chromosome prefix 'chr' is removed
      - only requested chromosomes are retained

    Returns:
      (gene_seg_df_normalized, gene_chunk_df_normalized_or_None, has_chunk_data_flag)
    """

    gdf = gene_seg_df.copy()
    if "chr" not in gdf.columns:
        raise ValueError("gene_seg_df must contain a 'chr' column.")
    gdf["chr"] = gdf["chr"].astype(str).str.replace("^chr", "", regex=True)
    gdf = gdf[gdf["chr"].isin(chr_order)].copy()

    gchunk = None
    if gene_chunk_df is not None and not gene_chunk_df.empty:
        gchunk = gene_chunk_df.copy()
        if "chr" not in gchunk.columns:
            raise ValueError("gene_chunk_df must contain a 'chr' column.")
        gchunk["chr"] = gchunk["chr"].astype(str).str.replace("^chr", "", regex=True)
        gchunk = gchunk[gchunk["chr"].isin(chr_order)].copy()

    has_chunk_data = gchunk is not None and not gchunk.empty
    return gdf, gchunk, has_chunk_data


def _normalize_targets_col(gdf: pd.DataFrame) -> tuple[pd.DataFrame, str | None]:
    """
    Normalize gene target count column naming.

    Accepts either:
      - 'n_targets'
      - 'n.targets'

    Returns:
      (normalized_dataframe, normalized_column_name_or_None)

    Does not modify values, only renames if needed.
    """

    out = gdf.copy()
    if "n_targets" in out.columns:
        return out, "n_targets"
    if "n.targets" in out.columns:
        out = out.rename(columns={"n.targets": "n_targets"})
        return out, "n_targets"
    return out, None


def _normalize_is_cancer_gene(gdf: pd.DataFrame) -> pd.DataFrame:
    """
    Normalize cancer gene annotation to boolean column.

    Creates:
      is_cancer_gene_bool

    Missing or null values default to False.
    """

    out = gdf.copy()
    if "is_cancer_gene" in out.columns:
        out["is_cancer_gene_bool"] = out["is_cancer_gene"].fillna(False).astype(bool)
    else:
        out["is_cancer_gene_bool"] = False
    return out


def _load_cnr_and_optional_pon(
    cnr_path: Path,
    chr_order: list[str],
    pon_path: Optional[Path],
    pct_spread: float,
    weight_thresh: float,
) -> tuple[pd.DataFrame, str, bool, float]:
    """
    Load CNR bins and optionally merge with PON spread information.

    Applies:
      - Chromosome normalization
      - Chromosome filtering
      - Optional PON merge (inner join behaviour)
      - Optional weight filtering
      - Global spread quantile threshold calculation

    Returns:
      merged_bins_df
      chromosome_column_name
      use_pon_flag
      spread_threshold_value

    Merged dataframe always contains a 'spread' column.
    """

    cnr = safe_read_csv(cnr_path, sep="\t")
    if cnr.empty:
        return pd.DataFrame(), "chromosome", False, float("inf")

    chr_col = "chromosome" if "chromosome" in cnr.columns else "chr"
    cnr[chr_col] = cnr[chr_col].astype(str).str.replace("^chr", "", regex=True)
    cnr = cnr[cnr[chr_col].isin(chr_order)].copy()

    required = {"start", "end", "log2"}
    missing = required - set(cnr.columns)
    if missing:
        raise ValueError(f"CNR missing required columns: {missing}")

    if "gene" not in cnr.columns:
        cnr["gene"] = ""

    use_pon = False
    merged = None

    if pon_path is not None and Path(pon_path).is_file():
        pon = safe_read_csv(pon_path, sep="\t")
        if not pon.empty and "spread" in pon.columns:
            try:
                merged = merge_cnr_pon(cnr, pon, chr_col=chr_col, spread_col="spread")
                if merged is not None and not merged.empty:
                    use_pon = True
            except Exception:
                use_pon = False

    if not use_pon:
        merged = cnr.copy()
        if "spread" not in merged.columns:
            merged["spread"] = np.nan

    merged = merged[merged[chr_col].isin(chr_order)].copy()

    # Apply global weight filter now if present (so spread quantile reflects used bins)
    if "weight" in merged.columns:
        merged = merged[merged["weight"] > weight_thresh].copy()

    spread_thresh = (
        float(merged["spread"].quantile(pct_spread)) if use_pon else float("inf")
    )
    return merged, chr_col, use_pon, spread_thresh


def _compute_row_flags(gdf: pd.DataFrame) -> pd.DataFrame:
    """
    Compute per-row CNV / LOH / PON significance flags.

    Adds:
      is_loh_or_cnv
      is_pon_signif
      cnv_flag = union of above

    Supports legacy and current column naming conventions.
    """

    out = gdf.copy()

    def _is_loh_or_cnv_row(row: pd.Series) -> bool:
        if "loh_flag" in row.index and pd.notna(row["loh_flag"]):
            if str(row["loh_flag"]).strip().upper() == "TRUE":
                return True
        for col in ("cnvkit_cnv_call", "purecn_cnv_call"):
            if col in row.index and pd.notna(row[col]):
                val = str(row[col]).strip().upper()
                if val in ("DELETION", "AMPLIFICATION"):
                    return True
        return False

    def _is_pon_signif_row(row: pd.Series) -> bool:
        if "pon_chunk_significance" in row.index and pd.notna(
            row["pon_chunk_significance"]
        ):
            return str(row["pon_chunk_significance"]).strip().lower() == "significant"

        if "pon_chunk_call" in row.index and pd.notna(row["pon_chunk_call"]):
            val = str(row["pon_chunk_call"]).strip().upper()
            if val in ("AMPLIFICATION", "DELETION"):
                return True

        if "pon_gene_call" in row.index and pd.notna(row["pon_gene_call"]):
            return str(row["pon_gene_call"]).strip().lower() == "significant"

        if "pon_gene_cnv_call" in row.index and pd.notna(row["pon_gene_cnv_call"]):
            val = str(row["pon_gene_cnv_call"]).strip().upper()
            if val in ("AMPLIFICATION", "DELETION"):
                return True

        if "pon_cnv_call" in row.index and pd.notna(row["pon_cnv_call"]):
            val = str(row["pon_cnv_call"]).strip().upper()
            if val in ("AMPLIFICATION", "DELETION"):
                return True

        if "pon_call" in row.index and pd.notna(row["pon_call"]):
            return str(row["pon_call"]).strip().lower() == "significant"

        return False

    out["is_loh_or_cnv"] = out.apply(_is_loh_or_cnv_row, axis=1)
    out["is_pon_signif"] = out.apply(_is_pon_signif_row, axis=1)
    out["cnv_flag"] = out["is_loh_or_cnv"] | out["is_pon_signif"]
    return out


def _compute_gene_level_highlights(
    gdf: pd.DataFrame,
    targets_col: str | None,
    highlight_only_cancer: bool,
    has_chunk_data: bool,
    min_gene_targets: int = 3,
    min_gene_targets_cancer: int = 5,
) -> pd.DataFrame | None:
    """
    Aggregate gene-level signal and determine highlighting eligibility.

    Computes per gene:
      - LOH/CNV presence
      - PON significance
      - Cancer gene flag
      - Total target count
      - highlight_gene decision

    Respects:
      - highlight_only_cancer mode
      - stricter target thresholds for cancer-only highlighting
      - fallback highlighting when chunk-level PON data is absent

    Returns gene-level summary dataframe or None if gene symbols unavailable.
    """

    if "gene.symbol" not in gdf.columns:
        return None

    grouped = gdf.groupby(["chr", "gene.symbol"], as_index=False)
    agg = {
        "has_loh_or_cnv": ("is_loh_or_cnv", "any"),
        "has_pon_sig": ("is_pon_signif", "any"),
        "is_cancer_gene": ("is_cancer_gene_bool", "any"),
    }
    if targets_col is not None:
        agg["total_targets"] = (targets_col, "sum")

    gene_level = grouped.agg(**agg)
    if "total_targets" not in gene_level.columns:
        gene_level["total_targets"] = 0.0

    gene_level["total_targets"] = gene_level["total_targets"].fillna(0.0)
    gene_level["is_cancer_gene"] = (
        gene_level["is_cancer_gene"].fillna(False).astype(bool)
    )

    min_req = min_gene_targets_cancer if highlight_only_cancer else min_gene_targets

    mask_base = (gene_level["has_loh_or_cnv"] | gene_level["has_pon_sig"]) & (
        gene_level["total_targets"] >= min_req
    )

    pon_only = gene_level["has_pon_sig"] & ~gene_level["has_loh_or_cnv"]
    mask_pon_cancer = ~pon_only | (pon_only & gene_level["is_cancer_gene"])

    mask_highlight = mask_base & mask_pon_cancer
    if highlight_only_cancer:
        mask_highlight &= gene_level["is_cancer_gene"]

    if not has_chunk_data:
        mask_highlight = mask_highlight | (
            gene_level["is_cancer_gene"] & (gene_level["total_targets"] >= min_req)
        )

    gene_level["highlight_gene"] = mask_highlight
    return gene_level


def _focus_genes_set(focus_genes: Optional[Collection[str]]) -> set[str] | None:
    """
    Normalize focus_genes collection into a string set.

    Returns None if focus_genes is None.
    """

    if focus_genes is None:
        return None
    return {str(g) for g in focus_genes}


def _restrict_to_focus_window(
    *,
    sub_bins: pd.DataFrame,
    g_chr: pd.DataFrame,
    focus_genes_set: set[str],
    focus_padding_bp: int,
) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """
    Restrict bins and gene rows to genomic window surrounding focus genes.

    Priority coordinate sources:
      1. region_start / region_end
      2. seg_start / seg_end
      3. fallback to CNR bin positions

    Returns:
      restricted_bins_df
      restricted_gene_rows_df
      list_of_focus_genes_present_on_chromosome
    """

    if g_chr.empty or "gene.symbol" not in g_chr.columns:
        return pd.DataFrame(), pd.DataFrame(), []

    g_focus_chr = g_chr[g_chr["gene.symbol"].isin(focus_genes_set)].copy()
    if g_focus_chr.empty:
        return pd.DataFrame(), pd.DataFrame(), []

    focus_genes_chr = sorted(set(g_focus_chr["gene.symbol"].astype(str).tolist()))

    # choose coordinate source
    if {"region_start", "region_end"}.issubset(g_focus_chr.columns):
        start_min = g_focus_chr["region_start"].min()
        end_max = g_focus_chr["region_end"].max()
    elif {"seg_start", "seg_end"}.issubset(g_focus_chr.columns):
        start_min = g_focus_chr["seg_start"].min()
        end_max = g_focus_chr["seg_end"].max()
    else:
        bins_focus = sub_bins[sub_bins["gene"].isin(focus_genes_chr)]
        if bins_focus.empty:
            return pd.DataFrame(), pd.DataFrame(), []
        start_min = bins_focus["start"].min()
        end_max = bins_focus["end"].max()

    start_min = max(0, int(start_min) - int(focus_padding_bp))
    end_max = int(end_max) + int(focus_padding_bp)

    sub_bins = sub_bins[
        (sub_bins["end"] >= start_min) & (sub_bins["start"] <= end_max)
    ].copy()
    if sub_bins.empty:
        return pd.DataFrame(), pd.DataFrame(), []

    if {"region_start", "region_end"}.issubset(g_chr.columns):
        g_chr = g_chr[
            (g_chr["region_end"] >= start_min) & (g_chr["region_start"] <= end_max)
        ].copy()
    elif {"seg_start", "seg_end"}.issubset(g_chr.columns):
        g_chr = g_chr[
            (g_chr["seg_end"] >= start_min) & (g_chr["seg_start"] <= end_max)
        ].copy()

    return sub_bins, g_chr, focus_genes_chr


def _compute_variable_x(
    sub: pd.DataFrame,
    highlighted_genes: np.ndarray,
    anti_factor: float,
    neutral_target_factor: float,
) -> pd.DataFrame:
    """
    Construct pseudo-position x coordinate using variable bin widths.

    Width rules:
      - Antitarget bins → anti_factor
      - Highlighted gene bins → 1.0
      - Other target bins → neutral_target_factor

    Adds:
      type
      bin_width
      x_coord
    """

    out = sub.copy()
    out["type"] = np.where(out["gene"] == "Antitarget", "Antitarget", "Target")

    def _bin_width(row: pd.Series) -> float:
        if row["type"] == "Antitarget":
            return anti_factor
        if row["gene"] in highlighted_genes:
            return 1.0
        return neutral_target_factor

    out["bin_width"] = out.apply(_bin_width, axis=1)
    out["x_coord"] = out["bin_width"].cumsum() - out["bin_width"] / 2
    return out


def _trim_pseudo_position_for_focus(
    sub: pd.DataFrame,
    focus_genes_chr: list[str],
    pseudo_pad: float,
) -> pd.DataFrame:
    """
    Restrict pseudo-position domain around focus genes.

    If multiple focus genes:
      Keep bins inside bounding window ± pseudo_pad.

    If single focus gene:
      Keep bins within pseudo distance threshold.

    Returns trimmed dataframe.
    """

    if not focus_genes_chr:
        return sub

    focus_bins = sub[sub["gene"].isin(focus_genes_chr)]
    if focus_bins.empty:
        return sub

    all_x = sub["x_coord"].to_numpy()
    if len(focus_genes_chr) >= 2:
        x_min = float(focus_bins["x_coord"].min()) - pseudo_pad
        x_max = float(focus_bins["x_coord"].max()) + pseudo_pad
        mask = (all_x >= x_min) & (all_x <= x_max)
    else:
        focus_x = focus_bins["x_coord"].to_numpy()
        dist = np.min(np.abs(all_x[:, None] - focus_x[None, :]), axis=1)
        mask = dist <= pseudo_pad

    return sub.loc[mask].copy()


def _add_smoothing(sub: pd.DataFrame, use_pon: bool, window: int) -> pd.DataFrame:
    """
    Add rolling median smoothing to log2 and optional PON metrics.

    Adds:
      log2_smooth
      spread_smooth (if PON available)
      pon_log2_smooth (if available, else zero baseline)

    Uses centered rolling window.
    """

    out = sub.sort_values("x_coord").copy()
    out["log2_smooth"] = out["log2"].rolling(window=window, center=True).median()
    if use_pon:
        out["spread_smooth"] = (
            out["spread"].rolling(window=window, center=True).median()
        )
        if "pon_log2" in out.columns:
            out["pon_log2_smooth"] = (
                out["pon_log2"].rolling(window=window, center=True).median()
            )
        else:
            out["pon_log2_smooth"] = 0.0
    else:
        out["spread_smooth"] = np.nan
        out["pon_log2_smooth"] = 0.0
    return out


def _pos_to_xcoord_fn(sub: pd.DataFrame) -> callable:
    """
    Build genomic position → pseudo-position mapping function.

    Uses nearest bin left-edge lookup via searchsorted.

    Returned function is safe for out-of-range positions.
    """

    bin_starts = sub["start"].to_numpy()
    x_coords = sub["x_coord"].to_numpy()

    def pos_to_xcoord(pos: int) -> float:
        idx = np.searchsorted(bin_starts, pos, side="right") - 1
        if idx < 0:
            idx = 0
        if idx >= len(x_coords):
            idx = len(x_coords) - 1
        return float(x_coords[idx])

    return pos_to_xcoord


def _collect_segments_for_chr(
    g_chr: pd.DataFrame, chr_name: str, y_clip: float
) -> pd.DataFrame:
    """
    Aggregate unique CNV/LOH segments from gene-level rows for plotting.

    Supports:
      - seg_log2
      - LOH mean values
      - CNVkit calls

    Adds clipped values for plotting stability.

    Returns per-segment dataframe.
    """

    if g_chr.empty or not {"seg_start", "seg_end"}.issubset(g_chr.columns):
        return pd.DataFrame()

    agg: dict[str, tuple[str, str]] = {}
    if "seg_log2" in g_chr.columns:
        agg["seg_log2"] = ("seg_log2", "first")
    if "seg_cn" in g_chr.columns:
        agg["seg_C"] = ("seg_cn", "first")
    if "loh_seg_mean" in g_chr.columns:
        agg["loh_seg_mean"] = ("loh_seg_mean", "first")
    if "loh_C" in g_chr.columns:
        agg["loh_C"] = ("loh_C", "first")
    if "cnvkit_cnv_call" in g_chr.columns:
        agg["cnvkit_cnv_call"] = ("cnvkit_cnv_call", "first")

    segs = (
        g_chr.dropna(subset=["seg_start", "seg_end"])
        .groupby(["chr", "seg_start", "seg_end"], as_index=False)
        .agg(**agg)
    )
    segs = segs[segs["chr"] == chr_name].sort_values("seg_start").copy()

    if "seg_log2" in segs.columns:
        segs["seg_log2_clipped"] = segs["seg_log2"].clip(-y_clip, y_clip)
    if "loh_seg_mean" in segs.columns:
        segs["loh_seg_mean_clipped"] = segs["loh_seg_mean"].clip(-y_clip, y_clip)

    return segs


def _stable_gene_color_fn():
    """
    Create deterministic gene → color mapping function.

    Uses tab20 colormap hashed by gene name to ensure stable
    color assignment across plots and runs.
    """

    cmap_all = colormaps.get_cmap("tab20")

    def _stable_gene_color(gname: str):
        h = abs(hash(str(gname)))
        return cmap_all(h % cmap_all.N)

    return _stable_gene_color


def _make_gene_colors(highlighted_genes: np.ndarray) -> dict[str, tuple]:
    """
    Create sequential colormap assignment for highlighted genes.

    Uses resampled tab20 colormap sized to number of genes.

    Returns gene → RGBA color mapping dictionary.
    """

    gene_to_color: dict[str, tuple] = {}
    if highlighted_genes.size == 0:
        return gene_to_color
    cmap = colormaps.get_cmap("tab20").resampled(max(int(highlighted_genes.size), 1))
    for i, gname in enumerate(highlighted_genes):
        gene_to_color[str(gname)] = cmap(i)
    return gene_to_color


def _draw_background_bins(
    ax, sub: pd.DataFrame, highlighted_genes: np.ndarray, y_col: str
):
    """
    Draw background bin scatter layer for non-highlighted bins.

    Separates:
      - Target bins
      - Antitarget bins

    Used to provide visual density context behind highlighted genes.
    """

    if highlighted_genes.size > 0:
        non_cnv = sub[~sub["gene"].isin(highlighted_genes)]
    else:
        non_cnv = sub

    bg_targets = non_cnv[non_cnv["type"] == "Target"]
    bg_antis = non_cnv[non_cnv["type"] == "Antitarget"]

    if not bg_antis.empty:
        ax.scatter(
            bg_antis["x_coord"],
            bg_antis[y_col],
            s=3,
            alpha=0.38,
            color="lightgrey",
            edgecolors="black",
            linewidths=0.2,
            label="Antitarget bins",
        )

    if not bg_targets.empty:
        ax.scatter(
            bg_targets["x_coord"],
            bg_targets[y_col],
            s=4,
            alpha=0.5,
            color="tab:blue",
            edgecolors="black",
            linewidths=0.25,
            label="Target bins (no highlighted CNV gene)",
        )


def _draw_highlighted_bins(
    ax,
    sub: pd.DataFrame,
    highlighted_genes: np.ndarray,
    gene_to_color: dict[str, tuple],
    y_col: str,
):
    """
    Draw highlighted gene bin scatter layer.

    Uses gene-specific colors and heavier styling than background bins.
    """

    for gene in highlighted_genes:
        gsub = sub[sub["gene"] == gene]
        if gsub.empty:
            continue
        color = gene_to_color.get(str(gene), "black")
        ax.scatter(
            gsub["x_coord"],
            gsub[y_col],
            s=8,
            alpha=0.9,
            color=color,
            edgecolors="black",
            linewidths=0.3,
        )


def _draw_density_bar(ax, sub: pd.DataFrame, y_min: float, y_max: float):
    """
    Draw target / antitarget density bar at top of plot.

    Represents bin type distribution across pseudo-position axis.
    """

    bar_height = 0.07 * (y_max - y_min)
    bar_bottom = y_max - bar_height

    mask_t = sub["type"] == "Target"
    if mask_t.any():
        ax.bar(
            sub.loc[mask_t, "x_coord"],
            bar_height,
            bottom=bar_bottom,
            width=sub.loc[mask_t, "bin_width"],
            align="center",
            color="tab:blue",
            alpha=0.6,
            linewidth=0,
            zorder=1,
        )

    mask_a = sub["type"] == "Antitarget"
    if mask_a.any():
        ax.bar(
            sub.loc[mask_a, "x_coord"],
            bar_height,
            bottom=bar_bottom,
            width=sub.loc[mask_a, "bin_width"],
            align="center",
            color="lightgrey",
            alpha=0.6,
            linewidth=0,
            zorder=1,
        )


def _draw_pon_band(ax, sub: pd.DataFrame, x: pd.Series, y_clip: float):
    """
    Draw PON noise band around PON smoothed log2 baseline.

    Band = PON mean ± PON spread.
    Clipped to plotting y-range.
    """

    if (
        "spread_smooth" not in sub.columns
        or "pon_log2_smooth" not in sub.columns
        or not sub["spread_smooth"].notna().any()
    ):
        return

    band_center = sub["pon_log2_smooth"].fillna(0.0)
    band_half = sub["spread_smooth"].fillna(0.0)
    band_bottom = (band_center - band_half).clip(-y_clip, y_clip)
    band_top = (band_center + band_half).clip(-y_clip, y_clip)

    ax.fill_between(
        x,
        band_bottom,
        band_top,
        alpha=0.4,
        step="mid",
        label="PON noise band",
    )


def _draw_segments(ax, segs_chr: pd.DataFrame):
    """
    Draw CNVkit / LOH segment horizontal lines.

    Segment color determined by CNV call:
      Amplification → red
      Deletion → blue
      Otherwise → black
    """

    if segs_chr.empty:
        return

    for _, srow in segs_chr.iterrows():
        xs = srow["x_start"]
        xe = srow["x_end"]

        if "seg_log2_clipped" in srow.index and pd.notna(srow["seg_log2_clipped"]):
            y_seg = srow["seg_log2_clipped"]
        elif "loh_seg_mean_clipped" in srow.index and pd.notna(
            srow["loh_seg_mean_clipped"]
        ):
            y_seg = srow["loh_seg_mean_clipped"]
        else:
            continue

        call_val = srow.get("cnvkit_cnv_call", np.nan)
        if pd.isna(call_val):
            seg_color = "black"
        else:
            c = str(call_val).strip().upper()
            if c == "AMPLIFICATION":
                seg_color = "red"
            elif c == "DELETION":
                seg_color = "royalblue"
            else:
                seg_color = "black"

        ax.hlines(y_seg, xs, xe, colors=seg_color, linewidth=1.8, alpha=0.9)


def _draw_chunk_pon_segments(
    ax,
    *,
    g_chunks_chr: pd.DataFrame,
    sub: pd.DataFrame,
    pos_to_xcoord: callable,
    chr_name: str,
    y_clip: float,
    highlight_only_cancer: bool,
    gene_level: pd.DataFrame | None,
    min_gene_targets_cancer: int,
) -> tuple[set[str], bool]:
    """
    Draw PON-driven CNV chunk segments as bright yellow lines.

    Also returns:
      - set of genes with PON CNV chunks
      - label_added flag (for legend control)

    Supports cancer-only filtering if gene-level summary provided.
    """

    pon_cnv_genes: set[str] = set()
    label_added = False

    if g_chunks_chr is None or g_chunks_chr.empty:
        return pon_cnv_genes, label_added

    if not {"region_start", "region_end", "mean_log2"}.issubset(g_chunks_chr.columns):
        return pon_cnv_genes, label_added

    start_span = sub["start"].min()
    end_span = sub["end"].max()
    chunks_span = g_chunks_chr[
        (g_chunks_chr["region_end"] >= start_span)
        & (g_chunks_chr["region_start"] <= end_span)
    ].copy()

    if "pon_chunk_call" not in chunks_span.columns:
        return pon_cnv_genes, label_added

    chunks_span["pon_chunk_call_norm"] = (
        chunks_span["pon_chunk_call"].astype(str).str.strip().str.upper()
    )
    chunks_span = chunks_span[
        chunks_span["pon_chunk_call_norm"].isin(["AMPLIFICATION", "DELETION"])
    ]

    if chunks_span.empty:
        return pon_cnv_genes, label_added

    for _, crow in chunks_span.iterrows():
        xs = pos_to_xcoord(int(crow["region_start"]))
        xe = pos_to_xcoord(int(crow["region_end"]))
        y = float(np.clip(float(crow["mean_log2"]), -y_clip, y_clip))

        line = ax.hlines(
            y,
            xs,
            xe,
            colors="yellow",
            linewidth=1.2,
            alpha=0.95,
            linestyles="solid",
            label="PON-based CNV segment" if not label_added else None,
        )
        line.set_path_effects(
            [pe.Stroke(linewidth=1.6, foreground="black"), pe.Normal()]
        )
        label_added = True

    if "gene.symbol" in chunks_span.columns:
        pon_cnv_genes = set(chunks_span["gene.symbol"].dropna().astype(str).tolist())

    if highlight_only_cancer and gene_level is not None and pon_cnv_genes:
        allowed = gene_level[
            (gene_level["chr"] == chr_name)
            & (gene_level["is_cancer_gene"])
            & (gene_level["total_targets"] >= min_gene_targets_cancer)
        ]["gene.symbol"].astype(str)
        pon_cnv_genes &= set(allowed)

    return pon_cnv_genes, label_added


def _draw_gene_labels(
    ax,
    *,
    g_chr: pd.DataFrame,
    sub: pd.DataFrame,
    label_genes: list[str],
    gene_to_color: dict[str, tuple],
    stable_color_fn: callable,
    pos_to_xcoord: callable,
    y_max: float,
    base_label_offset: float,
):
    """
    Draw gene boundary markers and rotated gene name labels.

    Determines genomic span using:
      region coordinates → preferred
      segment coordinates → fallback
      bin extents → final fallback
    """

    for gname in label_genes:
        color = gene_to_color.get(gname) or stable_color_fn(gname)

        g_rows = (
            g_chr[g_chr["gene.symbol"] == gname]
            if "gene.symbol" in g_chr.columns
            else pd.DataFrame()
        )
        if (not g_rows.empty) and {"region_start", "region_end"}.issubset(
            g_rows.columns
        ):
            g_start = int(g_rows["region_start"].min())
            g_end = int(g_rows["region_end"].max())
        elif (not g_rows.empty) and {"seg_start", "seg_end"}.issubset(g_rows.columns):
            g_start = int(g_rows["seg_start"].min())
            g_end = int(g_rows["seg_end"].max())
        else:
            bins_gene = sub[sub["gene"] == gname]
            if bins_gene.empty:
                continue
            g_start = int(bins_gene["start"].min())
            g_end = int(bins_gene["end"].max())

        xs = pos_to_xcoord(g_start)
        xe = pos_to_xcoord(g_end)

        for xg in (xs, xe):
            ax.axvline(xg, color=color, linewidth=1.0, alpha=0.95, zorder=4)

        x_mid = (xs + xe) / 2.0
        y_label = y_max - (base_label_offset * 0.5)
        ax.text(
            x_mid,
            y_label,
            gname,
            rotation=90,
            fontsize=9,
            ha="center",
            va="top",
            color=color,
            path_effects=[pe.Stroke(linewidth=1.5, foreground="white"), pe.Normal()],
        )


def _plot_baf_panel(ax, baf_chr: pd.DataFrame):
    """
    Draw BAF scatter panel and expected allele fraction reference lines.

    Includes:
      - 0.5 heterozygous reference
      - 1/3 and 2/3 optional imbalance references
    """

    if not baf_chr.empty:
        ax.scatter(baf_chr["x_coord"], baf_chr["BAF"], s=6, alpha=0.5)
    ax.axhline(0.5, color="gray", linewidth=0.8, linestyle="--")
    for frac in [1 / 3, 2 / 3]:
        ax.axhline(frac, color="lightgray", linewidth=0.6, linestyle=":")
    ax.set_ylim(0, 1)
    ax.set_ylabel("BAF")


def _output_png_path(outdir: Path, chr_name: str, focus_genes_chr: list[str]) -> Path:
    """
    Construct output PNG path for chromosome plot.

    Includes focus gene tag in filename if focus genes are active.
    """
    if focus_genes_chr:
        tag = "genes_" + "_".join(focus_genes_chr)
        return outdir / f"cnv_chr{chr_name}_{tag}_segments.png"
    return outdir / f"cnv_chr{chr_name}_segments.png"


# =============================================================================
# Main plotting function (now much thinner)
# =============================================================================


def plot_chromosomes(
    cnr_path: Path,
    vcf_path: Path,
    gene_seg_df: pd.DataFrame,
    outdir: Path,
    case_id: str,
    pon_path: Optional[Path] = None,
    gene_chunk_df: Optional[pd.DataFrame] = None,
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
      - PON spread band (if available)
      - smoothed log2 (CNR)
      - segments from gene_seg_df
      - PON CNV chunks from gene_chunk_df (if provided)
      - gene highlighting + labels
      - BAF from VCF
    """
    outdir.mkdir(parents=True, exist_ok=True)

    MIN_GENE_TARGETS = 3
    MIN_GENE_TARGETS_CANCER = 5
    PSEUDO_PAD = 50.0

    chr_order = _as_chr_order(include_y)

    # Normalize gene dfs
    gdf, gchunk, has_chunk_data = _normalize_gene_tables(
        gene_seg_df, chr_order, gene_chunk_df
    )
    gdf, targets_col = _normalize_targets_col(gdf)
    gdf = _normalize_is_cancer_gene(gdf)

    # Load VCF BAF
    vcf = load_vcf_with_baf(vcf_path=vcf_path, chr_order=chr_order)

    # Load CNR and optional PON
    merged, chr_col, use_pon, spread_thresh = _load_cnr_and_optional_pon(
        cnr_path=cnr_path,
        chr_order=chr_order,
        pon_path=pon_path,
        pct_spread=pct_spread,
        weight_thresh=weight_thresh,
    )
    if merged.empty:
        return

    # Spread filter (applied per-chr later as well); keep spread col always
    if use_pon:
        merged = merged[merged["spread"] <= spread_thresh].copy()

    # Compute gene flags & gene-level highlight decisions
    gdf = _compute_row_flags(gdf)
    gene_level = _compute_gene_level_highlights(
        gdf,
        targets_col=targets_col,
        highlight_only_cancer=highlight_only_cancer,
        has_chunk_data=has_chunk_data,
        min_gene_targets=MIN_GENE_TARGETS,
        min_gene_targets_cancer=MIN_GENE_TARGETS_CANCER,
    )

    focus_set = _focus_genes_set(focus_genes)
    stable_color = _stable_gene_color_fn()

    y_clip = float(y_abs_max)
    y_lim_chr = (-y_clip, y_clip)

    for chr_name in chr_order:
        sub = merged[merged[chr_col] == chr_name].copy()
        if sub.empty:
            continue

        sub = sub.sort_values("start", kind="stable")
        if "weight" in sub.columns:
            sub = sub[sub["weight"] > weight_thresh].copy()
            if sub.empty:
                continue

        g_chr = gdf[gdf["chr"] == chr_name].copy()
        g_chunks_chr = (
            gchunk[gchunk["chr"] == chr_name].copy()
            if (gchunk is not None and not gchunk.empty)
            else pd.DataFrame()
        )

        focus_genes_chr: list[str] = []
        if focus_set is not None:
            sub, g_chr, focus_genes_chr = _restrict_to_focus_window(
                sub_bins=sub,
                g_chr=g_chr,
                focus_genes_set=focus_set,
                focus_padding_bp=focus_padding_bp,
            )
            if sub.empty:
                continue

        # Determine highlighted genes
        if gene_level is not None and "gene.symbol" in g_chr.columns:
            genes_in_view = g_chr["gene.symbol"].dropna().astype(str).unique()
            gene_level_chr = gene_level[
                (gene_level["chr"] == chr_name)
                & (gene_level["gene.symbol"].astype(str).isin(genes_in_view))
            ]
            highlighted = (
                gene_level_chr[gene_level_chr["highlight_gene"]]["gene.symbol"]
                .astype(str)
                .unique()
            )
        else:
            if highlight_only_cancer and "is_cancer_gene_bool" in g_chr.columns:
                mask_gene = g_chr["cnv_flag"] & g_chr["is_cancer_gene_bool"]
            else:
                mask_gene = g_chr["cnv_flag"]
            highlighted = (
                g_chr.loc[mask_gene, "gene.symbol"].dropna().astype(str).unique()
                if "gene.symbol" in g_chr.columns
                else np.array([])
            )

        if focus_genes_chr:
            highlighted = np.union1d(
                highlighted, np.array(focus_genes_chr, dtype=object)
            )

        gene_to_color = _make_gene_colors(highlighted)

        # Build pseudo-x
        sub = _compute_variable_x(
            sub,
            highlighted_genes=highlighted,
            anti_factor=anti_factor,
            neutral_target_factor=neutral_target_factor,
        )
        sub = _trim_pseudo_position_for_focus(sub, focus_genes_chr, PSEUDO_PAD)
        if sub.empty:
            continue

        # Smooth + clip
        sub = _add_smoothing(sub, use_pon=use_pon, window=window)
        sub["log2_clipped"] = sub["log2"].clip(-y_clip, y_clip)
        sub["log2_smooth_clipped"] = sub["log2_smooth"].clip(-y_clip, y_clip)

        pos_to_xcoord = _pos_to_xcoord_fn(sub)

        # BAF chr (restricted to span if focus)
        baf_chr = vcf[vcf["CHROM"] == chr_name].sort_values("POS").copy()
        if focus_genes_chr and not baf_chr.empty:
            baf_chr = baf_chr[
                (baf_chr["POS"] >= sub["start"].min())
                & (baf_chr["POS"] <= sub["end"].max())
            ].copy()
        if not baf_chr.empty:
            baf_chr["x_coord"] = baf_chr["POS"].apply(lambda p: pos_to_xcoord(int(p)))

        # Segments from g_chr
        segs_chr = _collect_segments_for_chr(g_chr, chr_name, y_clip=y_clip)
        if not segs_chr.empty:
            segs_chr = segs_chr.copy()
            segs_chr["x_start"] = segs_chr["seg_start"].apply(
                lambda p: pos_to_xcoord(int(p))
            )
            segs_chr["x_end"] = segs_chr["seg_end"].apply(
                lambda p: pos_to_xcoord(int(p))
            )

        # --- figure ---
        fig, (ax1, ax2) = plt.subplots(
            2, 1, figsize=(14, 6), sharex=True, gridspec_kw={"height_ratios": [2, 1]}
        )

        x = sub["x_coord"]

        _draw_background_bins(ax1, sub, highlighted, y_col="log2_clipped")
        if highlighted.size > 0:
            _draw_highlighted_bins(
                ax1, sub, highlighted, gene_to_color, y_col="log2_clipped"
            )

        y_min, y_max = y_lim_chr
        _draw_density_bar(ax1, sub, y_min, y_max)

        if use_pon:
            _draw_pon_band(ax1, sub, x, y_clip=y_clip)

        ax1.plot(
            x,
            sub["log2_smooth_clipped"],
            linewidth=1.5,
            alpha=0.9,
            color="tab:green",
            label=f"log2 (median {window} bins)",
        )

        _draw_segments(ax1, segs_chr)

        pon_cnv_genes, _ = _draw_chunk_pon_segments(
            ax1,
            g_chunks_chr=g_chunks_chr,
            sub=sub,
            pos_to_xcoord=pos_to_xcoord,
            chr_name=chr_name,
            y_clip=y_clip,
            highlight_only_cancer=highlight_only_cancer,
            gene_level=gene_level,
            min_gene_targets_cancer=MIN_GENE_TARGETS_CANCER,
        )

        ax1.axhline(0, color="black", linewidth=0.8)
        ax1.set_ylim(*y_lim_chr)
        ax1.set_ylabel("log2 / PON band")

        title_suffix = "log2 vs PON spread" if use_pon else "log2 (no PON available)"
        title = f"Chr {chr_name} – {title_suffix}: {case_id}"
        if focus_genes_chr:
            title += f"  (genes: {', '.join(focus_genes_chr)})"
        ax1.set_title(title + "\n")
        ax1.legend(loc="upper right", fontsize=8)

        label_genes = sorted(set(map(str, highlighted)) | set(map(str, pon_cnv_genes)))
        if label_genes:
            _draw_gene_labels(
                ax1,
                g_chr=g_chr,
                sub=sub,
                label_genes=label_genes,
                gene_to_color=gene_to_color,
                stable_color_fn=stable_color,
                pos_to_xcoord=pos_to_xcoord,
                y_max=y_max,
                base_label_offset=base_label_offset,
            )

        _plot_baf_panel(ax2, baf_chr)
        ax2.set_xlabel(
            "Pseudo-position (highlighted genes expanded, other bins compressed)"
        )

        plt.tight_layout()
        out_png = _output_png_path(outdir, chr_name, focus_genes_chr)
        plt.savefig(out_png, dpi=150)
        plt.close(fig)
