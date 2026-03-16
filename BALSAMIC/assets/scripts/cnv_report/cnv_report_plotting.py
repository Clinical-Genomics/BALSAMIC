from __future__ import annotations

# Standard library
from pathlib import Path
from typing import Iterable

# Third-party
import hashlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colormaps
from matplotlib import patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from cnv_io import load_vcf_with_vaf


def _pos_to_xcoord_fn(bins: pd.DataFrame) -> callable:
    """
    Build genomic position → pseudo-position mapping function.

    Uses nearest bin left-edge lookup via searchsorted.

    Returned function is safe for out-of-range positions.
    """

    bin_starts = bins["start"].to_numpy()
    x_coords = bins["x_coord"].to_numpy()

    def pos_to_xcoord(pos: int) -> float:
        idx = np.searchsorted(bin_starts, pos, side="right") - 1
        if idx < 0:
            idx = 0
        if idx >= len(x_coords):
            idx = len(x_coords) - 1
        return float(x_coords[idx])

    return pos_to_xcoord


def _stable_gene_color_fn(cmap_name: str = "tab20"):
    """
    Deterministic gene → color mapping, stable across runs and machines.

    Uses sha256(gene_name) → integer → colormap index.
    """
    cmap = colormaps.get_cmap(cmap_name)

    def _stable_gene_color(gname: str):
        s = str(gname).encode("utf-8")
        h = hashlib.sha256(s).digest()
        idx = int.from_bytes(h[:4], byteorder="big") % cmap.N
        return cmap(idx)

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


def _draw_background_bins(ax, bins: pd.DataFrame, y_col: str):
    """
    Draw background scatter for NON-highlighted bins (bin-level).
    Requires columns: x_coord, type, and y_col.
    """
    bg_backbone = bins[bins["type"] == "Backbone"]
    bg_targets = bins[bins["type"] == "Target"]

    if not bg_backbone.empty:
        ax.scatter(
            bg_backbone["x_coord"],
            bg_backbone[y_col],
            s=6,
            alpha=0.5,
            color="black",
            edgecolors="none",
            label="Backbone bins",
        )

    if not bg_targets.empty:
        ax.scatter(
            bg_targets["x_coord"],
            bg_targets[y_col],
            s=7,
            alpha=0.6,
            color="blue",
            edgecolors="black",
            linewidths=0.3,
            label="Target bins (not highlighted)",
        )


def _draw_highlighted_bins(
    ax,
    *,
    bins: pd.DataFrame,
    highlighted_genes: np.ndarray,
    gene_to_color: dict[str, tuple],
    genes_col: str = "genes",
    y_col_bins: str = "log2_clipped",
) -> None:
    """
    Draw highlighted bins ONCE per unique bin (bins-only).

    Coloring rule:
      - If a bin contains exactly ONE highlighted gene → use that gene's color.
      - If a bin contains >1 highlighted genes → draw as neutral gray.
      - Bins with 0 highlighted genes are not drawn here.

    Requirements:
      - bins contains: x_coord, y_col_bins, genes_col (list[str]), is_highlight_bin
    """
    if bins is None or bins.empty:
        return
    if "is_highlight_bin" not in bins.columns:
        return

    highlight_set = set(map(str, highlighted_genes))
    if not highlight_set:
        return

    # Pre-filter to bins that contain at least one highlighted gene.
    # This avoids applying the per-bin gene inspection to all bins.
    highlighted_bins = bins[bins["is_highlight_bin"]].copy()
    if highlighted_bins.empty:
        return

    def highlighted_genes_in_bin(gene_list) -> list[str]:
        if not isinstance(gene_list, list) or not gene_list:
            return []
        # only genes that are in highlighted set
        return [gene for gene in gene_list if str(gene) in highlight_set]

    highlighted_bins["hi_genes"] = highlighted_bins[genes_col].apply(
        highlighted_genes_in_bin
    )

    highlighted_bins["hi_count"] = highlighted_bins["hi_genes"].apply(len)

    single = highlighted_bins[highlighted_bins["hi_count"] == 1].copy()
    multi = highlighted_bins[highlighted_bins["hi_count"] > 1].copy()

    # --- single highlighted gene: color by that gene ---
    if not single.empty:
        single["sole_gene"] = single["hi_genes"].apply(lambda lst: str(lst[0]))
        for gene, df_g in single.groupby("sole_gene", sort=False):
            color = gene_to_color.get(gene, "black")
            ax.scatter(
                df_g["x_coord"],
                df_g[y_col_bins],
                s=9,
                alpha=0.9,
                color=color,
                edgecolors="black",
                linewidths=0.6,
                label=None,
            )

    # --- multiple highlighted genes: neutral ---
    if not multi.empty:
        ax.scatter(
            multi["x_coord"],
            multi[y_col_bins],
            s=9,
            alpha=0.9,
            color="gray",
            edgecolors="black",
            linewidths=0.6,
            label=None,
        )


def _draw_pon_bars(ax, sub: pd.DataFrame, x: pd.Series, y_clip: float):
    """
    Draw PON noise bars around PON log2 baseline.

    Bar = PON mean ± PON spread.
    Clipped to plotting y-range.
    """

    band_center = sub["pon_log2"].fillna(0.0)
    band_spread = sub["pon_spread"].fillna(0.0)
    band_bottom = (band_center - band_spread).clip(-y_clip, y_clip)
    band_top = (band_center + band_spread).clip(-y_clip, y_clip)

    ax.fill_between(
        x,
        band_bottom,
        band_top,
        alpha=0.4,
        step="mid",
        label="PON target spread",
    )


def _draw_segments(
    ax,
    segs_chr: pd.DataFrame,
    *,
    caller_col: str = "caller",
    call_col: str = "cnv_call",
    y_col: str = "log2_clipped",
    x_start_col: str = "x_start",
    x_end_col: str = "x_end",
) -> None:
    """
    Draw CNV segments as horizontal lines for BOTH CNVkit and PureCN.

    Expected columns in segs_chr:
      - caller (CNVkit/PureCN)
      - cnv_call (AMPLIFICATION/DELETION/NEUTRAL)
      - log2_clipped
      - x_start, x_end

    Color by cnv_call:
      Amplification → red
      Deletion → royalblue
      Neutral → black

    Caller style:
      CNVkit → solid, thicker
      PureCN → dashed, thinner
    """
    if segs_chr is None or segs_chr.empty:
        return

    df = segs_chr.copy()

    # Standardise caller values
    df[caller_col] = df[caller_col].astype("string").fillna("").str.strip().str.upper()

    # Standardise call values
    calls = df[call_col].astype("string").fillna("").str.strip().str.upper()

    # Compute colors vectorized
    colors = np.where(
        calls.eq("AMPLIFICATION"),
        "red",
        np.where(calls.eq("DELETION"), "royalblue", "black"),
    )
    df["_seg_color"] = colors

    # Styles per caller
    style_map = {
        "CNVKIT": dict(linestyle="solid", linewidth=1.8, alpha=0.9, zorder=3),
        "PURECN": dict(linestyle=(0, (4, 2)), linewidth=1.2, alpha=0.9, zorder=3),
    }

    # Draw caller-by-caller so styles are consistent
    for caller, sub in df.groupby(caller_col, sort=False):
        if caller not in style_map:
            style = dict(linestyle="solid", linewidth=1.2, alpha=0.8, zorder=3)
        else:
            style = style_map[caller]

        for _, srow in sub.iterrows():
            xs = srow[x_start_col]
            xe = srow[x_end_col]
            y = srow[y_col]
            if pd.isna(xs) or pd.isna(xe) or pd.isna(y):
                continue

            ax.hlines(y, xs, xe, colors=srow["_seg_color"], **style)


def _draw_generegion_pon_segments(
    ax,
    *,
    generegions: pd.DataFrame,
    pos_to_xcoord: callable,
    y_clip: float,
    start_col: str = "region_start",
    end_col: str = "region_end",
    y_col: str = "mean_log2",
    pon_ind_col: str = "pon_region_indication",
    pon_sig_col: str = "pon_region_signal",
) -> None:
    """
    Draw PON-driven region segments if (and only if) PON columns exist.

    Draw condition:
      - pon_region_indication in {"GAIN","LOSS"} OR
      - pon_region_signal == "strong"  (if that column exists)

    Returns: None
    """
    empty_result: tuple[set[str], bool] = (set(), False)

    if generegions is None or generegions.empty:
        return empty_result

    # If there's no PON info at all, do nothing.
    has_ind = pon_ind_col in generegions.columns
    has_sig = pon_sig_col in generegions.columns
    if not (has_ind or has_sig):
        return empty_result

    regions = generegions.copy()

    regions = regions.dropna(subset=[start_col, end_col, y_col]).copy()
    if regions.empty:
        return empty_result

    # PON GAIN / LOSS signal filter (only strong signal)

    ind = regions[pon_ind_col].astype("string").str.upper().isin({"GAIN", "LOSS"})
    sig = regions[pon_sig_col].astype("string").str.lower().eq("strong")

    regions = regions[ind | sig].copy()

    if regions.empty:
        return empty_result

    # Draw
    label_text = "PON-based gain/loss indicator"
    label_added = False

    for _, region_row in regions.iterrows():
        xs = pos_to_xcoord(int(region_row[start_col]))
        xe = pos_to_xcoord(int(region_row[end_col]))
        y = float(np.clip(float(region_row[y_col]), -y_clip, y_clip))

        line = ax.hlines(
            y,
            xs,
            xe,
            colors="yellow",
            linewidth=1.2,
            alpha=0.95,
            linestyles="solid",
            label=label_text if not label_added else None,
        )
        line.set_path_effects(
            [pe.Stroke(linewidth=1.8, foreground="black"), pe.Normal()]
        )
        label_added = True


def draw_gene_labels_from_spans(
    ax,
    *,
    gene_spans: pd.DataFrame,
    label_genes: list[str],
    gene_to_color: dict[str, tuple],
    pos_to_xcoord: callable,
    y_max: float,
    base_label_offset: float,
    gene_col: str = "gene.symbol",
    start_col: str = "region_start",
    end_col: str = "region_end",
) -> None:
    if gene_spans is None or gene_spans.empty or not label_genes:
        return

    # Index spans by gene name for fast lookup
    gene_spans_indexed = gene_spans.copy()
    gene_spans_indexed[gene_col] = (
        gene_spans_indexed[gene_col].astype("string").str.strip()
    )
    gene_spans_indexed = gene_spans_indexed.set_index(gene_col, drop=True)

    for gene_name in label_genes:
        gene_name = str(gene_name).strip()

        if gene_name not in gene_spans_indexed.index:
            continue

        gene_start = int(gene_spans_indexed.at[gene_name, start_col])
        gene_end = int(gene_spans_indexed.at[gene_name, end_col])

        color = gene_to_color[gene_name]

        x_start = pos_to_xcoord(gene_start)
        x_end = pos_to_xcoord(gene_end)

        # draw vertical boundaries for gene span
        for x_boundary in (x_start, x_end):
            ax.axvline(x_boundary, color=color, linewidth=1.0, alpha=0.95, zorder=4)

        x_midpoint = (x_start + x_end) / 2.0
        y_label = y_max - (base_label_offset * 0.5)

        ax.text(
            x_midpoint,
            y_label,
            gene_name,
            rotation=90,
            fontsize=9,
            ha="center",
            va="top",
            color=color,
            path_effects=[pe.Stroke(linewidth=1.5, foreground="white"), pe.Normal()],
        )


def _plot_vaf_panel(ax, vaf_chr: pd.DataFrame):
    """
    Draw VAF scatter panel and expected allele fraction reference lines.

    Includes:
      - 0.5 heterozygous reference
      - 1/3 and 2/3 optional imbalance references
    """

    if not vaf_chr.empty:
        ax.scatter(vaf_chr["x_coord"], vaf_chr["VAF"], s=6, alpha=0.5)
    ax.axhline(0.5, color="gray", linewidth=0.8, linestyle="--")
    for frac in [1 / 3, 2 / 3]:
        ax.axhline(frac, color="lightgray", linewidth=0.6, linestyle=":")
    ax.set_ylim(0, 1)
    ax.set_ylabel("VAF")


def _compute_variable_width_bins(
    bins: pd.DataFrame,
    *,
    neutral_target_factor: float,
    backbone_factor: float,
    backbone_label: str = "backbone",
    genes_col: str = "genes",
) -> pd.DataFrame:
    """
    Compute variable-width pseudo genomic coordinates for plotting bins.

    This function converts genomic bins into a compressed plotting coordinate
    system where different bin types occupy different visual widths.

    The goal is to visually emphasize bins belonging to highlighted genes
    while compressing neutral and backbone bins, allowing important regions
    to appear larger in the plot without removing other bins.

    Width assignment:
        highlighted bins  → width = 1.0
        neutral targets   → width = neutral_target_factor
        backbone bins     → width = backbone_factor

    Using these widths, the function computes a cumulative pseudo-position
    (`x_coord`) representing the visual center of each bin. This coordinate
    system is later used to map genomic positions to the plotting axis.

    Requirements:
        - `genes_col` contains a list of genes overlapping each bin
        - `is_highlight_bin` column already exists

    Returns:
        DataFrame with added columns:
            is_backbone_bin : bool
            type            : {"Backbone","Target"}
            bin_width       : visual bin width
            x_coord         : pseudo-position used for plotting
    """

    out = bins.copy()

    # Determine whether a bin represents a backbone region
    # (backbone bins contain the special gene label "backbone")
    def is_backbone(glist) -> bool:
        return isinstance(glist, list) and (backbone_label in glist)

    # Mark backbone bins
    out["is_backbone_bin"] = out[genes_col].apply(is_backbone)

    # Classify bins for plotting (used later for coloring/background scatter)
    out["type"] = np.where(out["is_backbone_bin"], "Backbone", "Target")

    # ------------------------------------------------------------------
    # Assign visual width to each bin.
    #
    # Highlighted bins → full width (1.0) so important genes expand.
    # Neutral targets  → compressed using neutral_target_factor.
    # Backbone bins    → compressed further using backbone_factor.
    # ------------------------------------------------------------------
    width = np.full(len(out), float(neutral_target_factor), dtype=float)

    # Compress backbone bins
    width[out["is_backbone_bin"].to_numpy()] = float(backbone_factor)

    # Expand bins that contain highlighted gene
    width[out["is_highlight_bin"].to_numpy()] = 1.0

    # Store visual bin width
    out["bin_width"] = width

    # ------------------------------------------------------------------
    # Compute pseudo genomic x-coordinate for plotting.
    #
    # Bins are laid out sequentially according to their visual width.
    # The coordinate represents the center of each bin.
    # ------------------------------------------------------------------
    out["x_coord"] = np.cumsum(width) - (width / 2.0)
    return out


def _collapse_bins_to_unique(
    df: pd.DataFrame,
    *,
    key_cols: tuple[str, str, str] = ("chr", "start", "end"),
    gene_col: str = "gene.symbol",
    value_aggs: dict[str, str] | None = None,
    out_gene_col: str = "genes",
) -> pd.DataFrame:
    """
    Collapse exploded bin table to one row per (chr,start,end).

    - Aggregates gene symbols into a sorted unique list (column `out_gene_col`)
    - Aggregates requested value columns using `value_aggs` (default: "first")

    Assumption: for a given bin, numeric values like log2/depth are identical across expanded rows.
    """
    if value_aggs is None:
        value_aggs = {}

    keys = list(key_cols)

    # gene list per bin
    genes = (
        df.groupby(keys, as_index=False)[gene_col]
        .agg(lambda s: sorted(set(s.dropna().astype(str))))
        .rename(columns={gene_col: out_gene_col})
    )

    # numeric / other columns per bin
    if value_aggs:
        vals = df.groupby(keys, as_index=False).agg(
            **{k: (k, v) for k, v in value_aggs.items()}
        )
        out = vals.merge(genes, on=keys, how="left", validate="one_to_one")
    else:
        out = genes

    return out


def _add_cnv_call_color_legend(ax, segs_chr: pd.DataFrame | None) -> None:
    """
    Legend for CNV call colors (Amp/Del/Neutral), independent of caller.
    Only shows entries present on this chromosome.
    """
    if segs_chr is None or segs_chr.empty or "cnv_call" not in segs_chr.columns:
        return

    calls = segs_chr["cnv_call"].astype("string").fillna("").str.strip().str.upper()

    handles = []
    if calls.eq("AMPLIFICATION").any():
        handles.append(Line2D([0], [0], color="red", lw=2.0, label="Amp"))
    if calls.eq("DELETION").any():
        handles.append(Line2D([0], [0], color="royalblue", lw=2.0, label="Del"))
    if (~calls.isin(["AMPLIFICATION", "DELETION"]) | calls.eq("")).any():
        handles.append(Line2D([0], [0], color="black", lw=2.0, label="Neutral"))

    if not handles:
        return

    legend = ax.legend(handles=handles, loc="lower right", fontsize=7, frameon=True)
    ax.add_artist(legend)


def _annotate_highlight_bins(
    bins: pd.DataFrame,
    *,
    highlighted_genes: Iterable[str],
    genes_col: str = "genes",
    out_col: str = "is_highlight_bin",
) -> pd.DataFrame:
    """
    Add a boolean column indicating whether any gene in a bin
    is in the highlighted gene set.
    """
    highlighted = {str(gene) for gene in highlighted_genes}

    def contains_highlighted_gene(genes: list[str]) -> bool:
        if not genes:
            return False
        return any(gene in highlighted for gene in genes)

    out = bins.copy()
    out[out_col] = out[genes_col].apply(contains_highlighted_gene)
    return out


def _is_amp_del(series: pd.Series) -> pd.Series:
    """Return True for AMPLIFICATION or DELETION calls."""
    values = series.astype("string").fillna("").str.strip().str.upper()
    return values.isin({"AMPLIFICATION", "DELETION"})


def _has_loh_type(series: pd.Series) -> pd.Series:
    """Return True if the annotation contains LOH."""
    values = series.astype("string").fillna("").str.strip().str.upper()
    return values.str.contains("LOH", regex=False)


def _add_region_highlight_evidence(
    regions: pd.DataFrame,
    *,
    cnvkit_call_col: str,
    purecn_call_col: str,
    purecn_type_col: str,
) -> pd.DataFrame:
    """
    Add per-row evidence columns for CNV, LOH, and combined CNV/LOH.
    """
    out = regions.copy()

    has_cnv = pd.Series(False, index=out.index)
    has_cnv |= _is_amp_del(out[cnvkit_call_col])

    if purecn_call_col in out.columns:
        has_cnv |= _is_amp_del(out[purecn_call_col])

    has_loh = pd.Series(False, index=out.index)
    if purecn_type_col in out.columns:
        has_loh |= _has_loh_type(out[purecn_type_col])

    out["_has_cnv"] = has_cnv
    out["_has_loh"] = has_loh
    out["_has_cnv_or_loh"] = has_cnv | has_loh
    return out


def _summarize_gene_highlights(
    regions: pd.DataFrame,
    *,
    chr_col: str,
    gene_col: str,
    targets_col: str,
    cancer_col: str,
    highlight_only_cancer: bool,
    min_gene_targets: int,
    min_gene_targets_cancer: int,
) -> pd.DataFrame:
    """
    Collapse region-level evidence to one row per gene and compute highlight status.
    """
    gene_summary = (
        regions.groupby([chr_col, gene_col], as_index=False)
        .agg(
            total_targets=(targets_col, "sum"),
            is_cancer_gene=(cancer_col, "any"),
            has_cnv=("_has_cnv", "any"),
            has_loh=("_has_loh", "any"),
            has_cnv_or_loh=("_has_cnv_or_loh", "any"),
        )
        .copy()
    )

    cancer_gene_selected = gene_summary["is_cancer_gene"] & (
        gene_summary["total_targets"] >= float(min_gene_targets_cancer)
    )

    noncancer_gene_selected = (
        (~gene_summary["is_cancer_gene"])
        & gene_summary["has_cnv_or_loh"]
        & (gene_summary["total_targets"] >= float(min_gene_targets))
    )

    gene_summary["highlight_gene"] = (
        cancer_gene_selected
        if highlight_only_cancer
        else (cancer_gene_selected | noncancer_gene_selected)
    )

    return gene_summary


def compute_highlighted_genes_from_generegions(
    generegions_all_chromosomes_df: pd.DataFrame,
    *,
    chr_name: str,
    highlight_only_cancer: bool,
    min_gene_targets: int = 4,
    min_gene_targets_cancer: int = 4,
    chr_col: str = "chr",
    gene_col: str = "gene.symbol",
    targets_col: str = "n.targets",
    cancer_col: str = "is_cancer_gene",
    cnvkit_call_col: str = "cnvkit_cnv_call",
    purecn_call_col: str = "purecn_cnv_call",
    purecn_type_col: str = "purecn_type",
) -> tuple[np.ndarray, pd.DataFrame]:
    """
    Decide which genes should be highlighted on one chromosome.

    Highlight evidence:
      - CNV: CNVkit or PureCN call is AMPLIFICATION or DELETION
      - LOH: PureCN type contains LOH

    Rules:
      - Cancer genes are highlighted if they meet the cancer target threshold
      - Non-cancer genes are highlighted only if they have CNV/LOH evidence
        and meet the non-cancer target threshold
      - If highlight_only_cancer is True, only the cancer-gene rule is used

    Returns:
      highlighted_genes:
          Sorted numpy array of highlighted gene names
      gene_summary:
          Per-gene summary table with evidence columns and highlight flag
    """
    generegions = generegions_all_chromosomes_df[
        generegions_all_chromosomes_df[chr_col].astype("string") == str(chr_name)
    ].copy()

    generegions[gene_col] = generegions[gene_col].astype("string").str.strip()
    generegions = generegions[
        generegions[gene_col].notna() & generegions[gene_col].ne("")
    ].copy()
    generegions = generegions[generegions[gene_col].ne("backbone")].copy()

    if generegions.empty:
        return np.array([], dtype=object), pd.DataFrame()

    generegions = _add_region_highlight_evidence(
        generegions,
        cnvkit_call_col=cnvkit_call_col,
        purecn_call_col=purecn_call_col,
        purecn_type_col=purecn_type_col,
    )

    gene_summary = _summarize_gene_highlights(
        generegions,
        chr_col=chr_col,
        gene_col=gene_col,
        targets_col=targets_col,
        cancer_col=cancer_col,
        highlight_only_cancer=highlight_only_cancer,
        min_gene_targets=min_gene_targets,
        min_gene_targets_cancer=min_gene_targets_cancer,
    )

    highlighted_genes = gene_summary.loc[gene_summary["highlight_gene"], gene_col]
    highlighted_genes = highlighted_genes.astype(str).dropna().unique()
    highlighted_genes = np.array(
        sorted(set(highlighted_genes) - {"backbone"}),
        dtype=object,
    )

    return highlighted_genes, gene_summary


def compute_gene_spans_from_generegions(
    generegions: pd.DataFrame,
    *,
    gene_col: str = "gene.symbol",
    start_col: str = "region_start",
    end_col: str = "region_end",
) -> pd.DataFrame:
    """
    Return one row per gene with full-span coordinates:
      region_start = min(region starts)
      region_end   = max(region ends)
    """
    if generegions is None or generegions.empty:
        return pd.DataFrame(columns=[gene_col, start_col, end_col])

    spans = generegions[[gene_col, start_col, end_col]].copy()
    spans[gene_col] = spans[gene_col].astype("string").str.strip()

    spans[start_col] = pd.to_numeric(spans[start_col], errors="coerce")
    spans[end_col] = pd.to_numeric(spans[end_col], errors="coerce")
    spans = spans.dropna(subset=[gene_col, start_col, end_col]).copy()

    gene_span = spans.groupby(gene_col, as_index=False).agg(
        **{start_col: (start_col, "min"), end_col: (end_col, "max")}
    )

    # ensure start <= end
    gene_span = gene_span[gene_span[end_col] >= gene_span[start_col]].copy()
    return gene_span


def get_loh_segments_for_chr(
    segments_df: pd.DataFrame,
    *,
    chr_name: str,
    chr_col: str = "chr",
    caller_col: str = "caller",
    start_col: str = "start",
    end_col: str = "end",
    purecn_type_col: str = "purecn_type",
) -> pd.DataFrame:
    df = segments_df.copy()

    df = df[df[chr_col].astype("string") == str(chr_name)]
    if df.empty or purecn_type_col not in df.columns or caller_col not in df.columns:
        return df.iloc[0:0].copy()

    is_purecn = df[caller_col].astype("string").str.upper().eq("PURECN")

    t = df[purecn_type_col].astype("string").fillna("").str.strip().str.upper()
    is_loh = t.str.contains("LOH", regex=False)

    # exclude unreliable LOH segments (determined by M.flagged)
    is_reliable = ~t.str.contains("UNRELIABLE", regex=False)

    out = df.loc[is_purecn & is_loh & is_reliable, [start_col, end_col]].copy()
    out[start_col] = pd.to_numeric(out[start_col], errors="coerce")
    out[end_col] = pd.to_numeric(out[end_col], errors="coerce")
    out = out.dropna(subset=[start_col, end_col])
    out = out.sort_values(start_col, kind="stable")
    return out


def draw_loh_spans_on_vaf(
    ax,
    loh_segs: pd.DataFrame,
    *,
    pos_to_xcoord: callable,
    start_col: str = "start",
    end_col: str = "end",
) -> bool:
    """
    Draw LOH intervals as translucent vertical shading on the VAF panel.

    Returns True if at least one interval was drawn.
    """
    if loh_segs is None or loh_segs.empty:
        return False

    drew = False
    for _, r in loh_segs.iterrows():
        xs = pos_to_xcoord(int(r[start_col]))
        xe = pos_to_xcoord(int(r[end_col]))
        if xe < xs:
            xs, xe = xe, xs

        ax.axvspan(xs, xe, alpha=0.12)
        drew = True

    return drew


def _add_segment_caller_style_legend(ax) -> None:
    handles = [
        Line2D([0], [0], color="black", lw=1.8, linestyle="solid", label="CNVkit"),
        Line2D([0], [0], color="black", lw=1.2, linestyle=(0, (4, 2)), label="PureCN"),
    ]
    leg = ax.legend(handles=handles, loc="lower left", fontsize=7, frameon=True)
    ax.add_artist(leg)


def plot_chromosomes(
    cnr_df: pd.DataFrame,
    vcf_path: Path,
    segments_df: pd.DataFrame,
    generegions_all_chromosomes_df: pd.DataFrame,
    outdir: Path,
    case_id: str,
    pon_df: pd.DataFrame | None = None,
    backbone_factor: float = 0.4,
    neutral_target_factor: float = 0.4,
    highlight_only_cancer: bool = False,
    window: int = 5,
    base_label_offset: float = 1.5,
    y_abs_max: float = 3.0,
) -> None:
    """
    Create one PNG per chromosome with:
      - PON spread band (if available)
      - smoothed log2 (CNR)
      - segments from gene_seg_df
      - PON CNV regions from generegion_df (if provided)
      - gene highlighting + labels
      - BAF from VCF
    """
    outdir.mkdir(parents=True, exist_ok=True)

    MIN_GENE_TARGETS = 4
    MIN_GENE_TARGETS_CANCER = 4

    KEY = ("chr", "start", "end")

    # --- collapse CNR to unique bins with gene list ---
    cnr_bins = _collapse_bins_to_unique(
        cnr_df,
        key_cols=KEY,
        value_aggs={
            "log2": "first",
            "depth": "first",
        },
        out_gene_col="genes",
    )

    # --- collapse PON to unique bins (no need for genes for PON, but harmless if present) ---
    if pon_df is not None:
        pon_bins = _collapse_bins_to_unique(
            pon_df,
            key_cols=KEY,
            value_aggs={
                "pon_log2": "first",
                "pon_spread": "first",
            },
            out_gene_col="pon_genes",  # keep separate if you want; you can also omit entirely
        )

        merged = cnr_bins.merge(
            pon_bins[list(KEY) + ["pon_log2", "pon_spread"]],
            on=list(KEY),
            how="left",
            validate="one_to_one",
        )
        use_pon = True
    else:
        merged = cnr_bins.copy()
        merged["pon_log2"] = np.nan
        merged["pon_spread"] = np.nan
        use_pon = False

    # Load VCF
    chr_order = [str(i) for i in range(1, 23)] + ["X", "Y"]
    vcf = load_vcf_with_vaf(vcf_path=vcf_path, chr_order=chr_order)

    merged = merged[merged["chr"].isin(chr_order)].copy()

    y_clip = float(y_abs_max)
    y_lim_chr = (-y_clip, y_clip)

    for chr_name in chr_order:
        highlighted, gene_summary = compute_highlighted_genes_from_generegions(
            generegions_all_chromosomes_df,
            chr_name=chr_name,
            highlight_only_cancer=highlight_only_cancer,
            min_gene_targets=MIN_GENE_TARGETS,
            min_gene_targets_cancer=MIN_GENE_TARGETS_CANCER,
        )

        # collapsed bin-level table for this chromosome
        chr_bins = merged[merged["chr"] == chr_name].copy()
        if chr_bins.empty:
            continue
        chr_bins = chr_bins.sort_values("start", kind="stable")
        if chr_bins.empty:
            continue

        generegions: pd.DataFrame = generegions_all_chromosomes_df[
            generegions_all_chromosomes_df["chr"] == chr_name
        ].copy()

        gene_to_color = _make_gene_colors(highlighted)

        # Mark which collapsed bins contain any highlighted gene
        chr_bins = _annotate_highlight_bins(
            chr_bins,
            highlighted_genes=highlighted,
            genes_col="genes",
            out_col="is_highlight_bin",
        )

        # Compute variable x-coordinates on collapsed bins
        x_coordinate_bins = _compute_variable_width_bins(
            chr_bins,
            neutral_target_factor=neutral_target_factor,
            backbone_factor=backbone_factor,
            backbone_label="backbone",
            genes_col="genes",
        )

        # Smooth + clip on UNIQUE bins
        x_coordinate_bins = x_coordinate_bins.sort_values(
            "x_coord", kind="stable"
        ).copy()
        x_coordinate_bins["log2_smooth"] = (
            x_coordinate_bins["log2"].rolling(window=window, center=True).median()
        )
        x_coordinate_bins["log2_clipped"] = x_coordinate_bins["log2"].clip(
            -y_clip, y_clip
        )
        x_coordinate_bins["log2_smooth_clipped"] = x_coordinate_bins[
            "log2_smooth"
        ].clip(-y_clip, y_clip)

        # Map function should use UNIQUE bins (stable mapping)
        pos_to_xcoord = _pos_to_xcoord_fn(x_coordinate_bins)

        # VAF chr
        vaf_chr = vcf[vcf["CHROM"] == chr_name].sort_values("POS").copy()
        if not vaf_chr.empty:
            vaf_chr["x_coord"] = vaf_chr["POS"].apply(lambda p: pos_to_xcoord(int(p)))

        # Segments from g_chr
        segs_chr = segments_df[
            segments_df["chr"].astype("string") == str(chr_name)
        ].copy()
        segs_chr["log2_clipped"] = segs_chr["log2"].clip(-y_clip, y_clip)
        if not segs_chr.empty:
            segs_chr = segs_chr.copy()
            segs_chr["x_start"] = segs_chr["start"].apply(
                lambda p: pos_to_xcoord(int(p))
            )
            segs_chr["x_end"] = segs_chr["end"].apply(lambda p: pos_to_xcoord(int(p)))

        # --- figure ---
        fig, (ax1, ax2) = plt.subplots(
            2, 1, figsize=(14, 6), sharex=True, gridspec_kw={"height_ratios": [2, 1]}
        )

        x = x_coordinate_bins["x_coord"]

        # background should be BIN-LEVEL and exclude bins that have any highlighted gene
        bg_bins = x_coordinate_bins[~x_coordinate_bins["is_highlight_bin"]].copy()
        _draw_background_bins(ax1, bg_bins, y_col="log2_clipped")

        # highlighted bins: one point per BIN, gene-colored if unique highlighted gene else neutral
        if highlighted.size > 0:
            _draw_highlighted_bins(
                ax1,
                bins=x_coordinate_bins,
                highlighted_genes=highlighted,
                gene_to_color=gene_to_color,
                genes_col="genes",
                y_col_bins="log2_clipped",
            )

        y_min, y_max = y_lim_chr

        if use_pon:
            _draw_pon_bars(ax1, x_coordinate_bins, x, y_clip=y_clip)

        ax1.plot(
            x,
            x_coordinate_bins["log2_smooth_clipped"],
            linewidth=1.5,
            alpha=0.9,
            color="tab:green",
            label=f"log2 (median {window} bins)",
        )

        _draw_generegion_pon_segments(
            ax1,
            generegions=generegions,
            pos_to_xcoord=pos_to_xcoord,
            y_clip=y_clip,
        )

        ax1.axhline(0, color="black", linewidth=0.8)
        ax1.set_ylim(*y_lim_chr)
        ax1.set_ylabel("log2 / PON band")

        title_suffix = "log2 vs PON spread" if use_pon else "log2 (no PON available)"
        title = f"Chr {chr_name} – {title_suffix}: {case_id}"
        ax1.set_title(title + "\n")
        main_leg = ax1.legend(loc="upper right", fontsize=8)

        gene_spans = compute_gene_spans_from_generegions(generegions)
        label_genes = set(map(str, highlighted))
        label_genes.discard("backbone")
        label_genes = sorted(label_genes)
        if label_genes:
            draw_gene_labels_from_spans(
                ax1,
                gene_spans=gene_spans,
                label_genes=label_genes,
                gene_to_color=gene_to_color,
                pos_to_xcoord=pos_to_xcoord,
                y_max=y_max,
                base_label_offset=base_label_offset,
            )

        _draw_segments(ax1, segs_chr)
        _add_segment_caller_style_legend(ax1)  # dashed vs solid
        _add_cnv_call_color_legend(ax1, segs_chr)  # red/blue/black

        ax1.add_artist(main_leg)
        loh_chr = get_loh_segments_for_chr(segments_df, chr_name=chr_name)
        # draw shading behind points:
        drew_loh = draw_loh_spans_on_vaf(ax2, loh_chr, pos_to_xcoord=pos_to_xcoord)

        _plot_vaf_panel(ax2, vaf_chr)

        # Add LOH legend handle if any LOH was drawn
        if drew_loh:
            loh_handle = Patch(
                alpha=0.12,
                label="PureCN LOH segment (reduced allelic balance expected)",
            )
            handles, labels = ax2.get_legend_handles_labels()
            ax2.legend(
                handles + [loh_handle],
                labels + [loh_handle.get_label()],
                loc="upper right",
                fontsize=7,
            )

        ax2.set_xlabel(
            "Pseudo-position (highlighted genes expanded, other bins compressed)"
        )

        plt.tight_layout()
        out_png = outdir / f"cnv_chr{chr_name}_segments.png"
        plt.savefig(out_png, dpi=150)
        plt.close(fig)
