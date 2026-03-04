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


from cnv_io import load_vcf_with_vaf


def _compute_row_flags(gdf: pd.DataFrame) -> pd.DataFrame:
    """
    Compute per-row CNV / LOH / PON significance flags.

    Logic:
      - is_loh_or_cnv:
            True if:
                * purecn 'type' column is non-empty
                OR
                * cnvkit_cnv_call in {"DELETION", "AMPLIFICATION"}
                OR
                * purecn_cnv_call in {"DELETION", "AMPLIFICATION"}

      - is_pon_signif:
            True if pon_chunk_significance == "significant"
            OR pon_chunk_indication in {"GAIN", "LOSS"}

      - cnv_flag = union of the above
    """

    out = gdf.copy()

    # ---------------------------------------------------------
    # LOH / CNV logic
    # ---------------------------------------------------------

    is_loh_or_cnv = pd.Series(False, index=out.index)

    # PureCN type present (non-empty)
    if "loh_flag" in out.columns:
        type_series = out["loh_flag"].astype("string").fillna("").str.strip()
        is_loh_or_cnv |= type_series.ne("")

    # CNVkit or PureCN explicit CNV calls
    for col in ("cnvkit_cnv_call", "purecn_cnv_call"):
        if col in out.columns:
            call_series = out[col].astype("string").fillna("").str.strip().str.upper()
            is_loh_or_cnv |= call_series.isin({"DELETION", "AMPLIFICATION"})

    out["is_loh_or_cnv"] = is_loh_or_cnv

    # ---------------------------------------------------------
    # PON significance logic
    # ---------------------------------------------------------

    is_pon_signif = pd.Series(False, index=out.index)

    if "pon_chunk_significance" in out.columns:
        signif_series = (
            out["pon_chunk_significance"]
            .astype("string")
            .fillna("")
            .str.strip()
            .str.lower()
        )
        is_pon_signif |= signif_series.eq("significant")

    if "pon_chunk_indication" in out.columns:
        indication_series = (
            out["pon_chunk_indication"]
            .astype("string")
            .fillna("")
            .str.strip()
            .str.upper()
        )
        is_pon_signif |= indication_series.isin({"GAIN", "LOSS"})

    out["is_pon_signif"] = is_pon_signif

    # ---------------------------------------------------------
    # Combined flag
    # ---------------------------------------------------------

    out["cnv_flag"] = out["is_loh_or_cnv"] | out["is_pon_signif"]

    return out


def _compute_gene_level_highlights(
    gdf: pd.DataFrame,
    targets_col: str,
    highlight_only_cancer: bool,
    min_gene_targets: int = 5,
    min_gene_targets_cancer: int = 5,
) -> pd.DataFrame:
    """
    Decide which genes should be highlighted.

    Rules:
      1) Cancer genes are ALWAYS highlighted if they meet the cancer target threshold.
      2) If highlight_only_cancer is False, also highlight non-cancer genes with LOH/CNV
         if they meet the non-cancer target threshold.
    """

    grouped = gdf.groupby(["chr", "gene.symbol"], as_index=False)
    gene_level = grouped.agg(
        has_loh_or_cnv=("is_loh_or_cnv", "any"),
        is_cancer_gene=("is_cancer_gene", "any"),
        total_targets=(targets_col, "sum"),
    )

    gene_level["total_targets"] = gene_level["total_targets"].fillna(0.0)
    gene_level["is_cancer_gene"] = (
        gene_level["is_cancer_gene"].fillna(False).astype(bool)
    )

    # 1) Cancer genes: always highlight if enough targets
    cancer_ok = gene_level["is_cancer_gene"] & (
        gene_level["total_targets"] >= min_gene_targets_cancer
    )

    # 2) Non-cancer CNV/LOH genes: optionally highlight if enough targets
    noncancer_cnv_ok = (
        (~gene_level["is_cancer_gene"])
        & gene_level["has_loh_or_cnv"]
        & (gene_level["total_targets"] >= min_gene_targets)
    )

    if highlight_only_cancer:
        mask_highlight = cancer_ok
    else:
        mask_highlight = cancer_ok | noncancer_cnv_ok

    gene_level["highlight_gene"] = mask_highlight
    return gene_level


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
    Aggregate unique CNV segments from gene-level rows for plotting.

    Supports:
      - seg_log2
      - CNVkit calls

    Adds clipped values for plotting stability.

    Returns per-segment dataframe.
    """

    agg: dict[str, tuple[str, str]] = {}

    agg["cnvkit_seg_raw_log2"] = ("cnvkit_seg_raw_log2", "first")
    agg["cnvkit_cnv_call"] = ("cnvkit_cnv_call", "first")

    segs = (
        g_chr.dropna(subset=["cnvkit_seg_start", "cnvkit_seg_end"])
        .groupby(["chr", "cnvkit_seg_start", "cnvkit_seg_end"], as_index=False)
        .agg(**agg)
    )
    segs = segs[segs["chr"] == chr_name].sort_values("cnvkit_seg_start").copy()

    segs["cnvkit_seg_log2_clipped"] = segs["cnvkit_seg_raw_log2"].clip(-y_clip, y_clip)

    return segs


def _stable_gene_color_fn(cmap_name: str = "tab20"):
    """
    Deterministic gene → color mapping, stable across runs and machines.

    Uses md5(gene_name) → integer → colormap index.
    """
    cmap = colormaps.get_cmap(cmap_name)

    def _stable_gene_color(gname: str):
        s = str(gname).encode("utf-8")
        h = hashlib.md5(s).hexdigest()  # 32 hex chars
        idx = int(h[:8], 16) % cmap.N  # use first 32 bits
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

    hi_set = set(map(str, highlighted_genes))
    if not hi_set:
        return

    # Only consider bins that are flagged as highlighted (fast path)
    b = bins[bins["is_highlight_bin"]].copy()
    if b.empty:
        return

    def hi_genes_in_bin(glist) -> list[str]:
        if not isinstance(glist, list) or not glist:
            return []
        # only genes that are in highlighted set
        return [g for g in glist if str(g) in hi_set]

    b["hi_genes"] = b[genes_col].apply(hi_genes_in_bin)

    single = b[b["hi_genes"].apply(len) == 1].copy()
    multi = b[b["hi_genes"].apply(len) > 1].copy()

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
    Draw PON noise bars around PON smoothed log2 baseline.

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

        if "cnvkit_seg_log2_clipped" in srow.index and pd.notna(
            srow["cnvkit_seg_log2_clipped"]
        ):
            y_seg = srow["cnvkit_seg_log2_clipped"]
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
    gene_level: pd.DataFrame,
    min_gene_targets_cancer: int,
) -> tuple[set[str], bool]:
    """
    Draw PON-driven CNV chunk segments (gain/loss) as yellow horizontal lines.

    Returns:
      (pon_cnv_genes, label_added)

    Notes:
      - A chunk is drawn if it overlaps the plotted genomic span and has pon_chunk_indication in {"GAIN","LOSS"}.
      - If highlight_only_cancer=True, returned genes are restricted to cancer genes with enough targets
        (drawing is unchanged; this only affects which genes may be labeled later).
    """
    empty_result: tuple[set[str], bool] = (set(), False)

    if g_chunks_chr is None or g_chunks_chr.empty:
        return empty_result
    if sub.empty:
        return empty_result

    # 1) Keep only chunks that overlap the span of the bins we are plotting
    span_start = int(sub["start"].min())
    span_end = int(sub["end"].max())

    chunks = g_chunks_chr.loc[
        (g_chunks_chr["region_end"] >= span_start)
        & (g_chunks_chr["region_start"] <= span_end)
    ].copy()

    if chunks.empty:
        return empty_result

    # 2) Keep only chunks with a PON gain/loss call
    pon_col = "pon_chunk_indication"
    if pon_col not in chunks.columns:
        return empty_result

    chunks[pon_col] = chunks[pon_col].astype(str).str.strip().str.upper()
    chunks = chunks[chunks[pon_col].isin({"GAIN", "LOSS"})]

    if chunks.empty:
        return empty_result

    # 3) Draw all matching chunks
    label_text = "PON-based gain/loss indicator"
    label_added = False

    for _, crow in chunks.iterrows():
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
            label=label_text if not label_added else None,
        )
        line.set_path_effects(
            [pe.Stroke(linewidth=1.8, foreground="black"), pe.Normal()]
        )
        label_added = True

    # 4) Collect genes represented by these chunks (used later for labeling)
    pon_cnv_genes: set[str] = set()
    if "gene.symbol" in chunks.columns:
        pon_cnv_genes = set(chunks["gene.symbol"].dropna().astype(str).tolist())

    # 5) Optionally restrict returned gene set (does NOT affect drawing)
    if highlight_only_cancer:
        if pon_cnv_genes:
            allowed = (
                gene_level.loc[
                    (gene_level["chr"] == chr_name)
                    & (gene_level["is_cancer_gene"])
                    & (gene_level["total_targets"] >= min_gene_targets_cancer),
                    "gene.symbol",
                ]
                .dropna()
                .astype(str)
            )

            pon_cnv_genes &= set(allowed.tolist())

    return pon_cnv_genes, label_added


def _draw_gene_labels(
    ax,
    *,
    g_chr: pd.DataFrame,
    label_genes: list[str],
    gene_to_color: dict[str, tuple],
    stable_color_fn: callable,
    pos_to_xcoord: callable,
    y_max: float,
    base_label_offset: float,
) -> None:
    """
    Draw gene boundary markers and rotated gene name labels.

    Requires `g_chr` to contain, per gene:
      - gene.symbol
      - region_start
      - region_end

    If a gene is missing coordinates, it is skipped.
    """
    if g_chr is None or g_chr.empty:
        return

    spans = g_chr.loc[
        g_chr["gene.symbol"].notna(), ["gene.symbol", "region_start", "region_end"]
    ].copy()
    # Ensure numeric-ish (won't crash if already int)
    spans["region_start"] = pd.to_numeric(spans["region_start"], errors="coerce")
    spans["region_end"] = pd.to_numeric(spans["region_end"], errors="coerce")

    gene_span = (
        spans.dropna(subset=["region_start", "region_end"])
        .groupby("gene.symbol", as_index=True)
        .agg(region_start=("region_start", "min"), region_end=("region_end", "max"))
    )

    for gname in label_genes:
        if gname not in gene_span.index:
            continue

        g_start = int(gene_span.at[gname, "region_start"])
        g_end = int(gene_span.at[gname, "region_end"])
        if g_end < g_start:
            continue

        color = gene_to_color.get(gname) or stable_color_fn(gname)

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


def _compute_variable_x_binlevel_OLD(
    sub: pd.DataFrame,
    highlighted_genes: np.ndarray,
    neutral_target_factor: float,
    backbone_factor: float = 0.3,
    backbone_label: str = "backbone",
    key_cols: tuple[str, str, str] = ("chr", "start", "end"),
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Compute pseudo-x on UNIQUE bins (chr,start,end) so exploded multi-gene rows
    do not inflate width.

    Returns:
      sub_full: original sub with x_coord/bin_width/type/is_highlight_bin merged in
      bins:     unique-bin dataframe (one row per bin) with x_coord/bin_width/etc
               (use this for smoothing/lines/background bins)

    Requirements in `sub`:
      - key_cols (default chr/start/end)
      - gene.symbol
      - log2
    Optional in `sub`:
      - pon_log2
      - pon_spread
    """
    out = sub.copy()

    # Ensure gene is a string dtype (may contain <NA>)
    gene = out["gene.symbol"].astype("string")

    # Set of highlighted gene names (string)
    hi = set(map(str, highlighted_genes))

    # Add per-row helper flags
    tmp = out.assign(
        _is_backbone=gene.eq(backbone_label),
        _is_hi=gene.fillna("").astype(str).isin(hi),
    )

    # Build aggregation dictionary safely (only include PON fields if present)
    agg: dict[str, tuple[str, str]] = {
        "is_backbone_bin": ("_is_backbone", "any"),
        "is_highlight_bin": ("_is_hi", "any"),
        "log2": ("log2", "first"),
    }
    if "pon_log2" in tmp.columns:
        agg["pon_log2"] = ("pon_log2", "first")
    if "pon_spread" in tmp.columns:
        agg["pon_spread"] = ("pon_spread", "first")

    # Aggregate to one row per unique bin
    g = tmp.groupby(list(key_cols), as_index=False).agg(**agg)

    # Bin-level plotting type
    g["type"] = np.where(g["is_backbone_bin"], "Backbone", "Target")

    # Bin-level width: neutral by default, compress backbone, expand highlighted bins
    width = np.full(len(g), float(neutral_target_factor), dtype=float)
    width[g["is_backbone_bin"].to_numpy()] = float(backbone_factor)
    width[g["is_highlight_bin"].to_numpy()] = 1.0

    g["bin_width"] = width
    g["x_coord"] = np.cumsum(width) - (width / 2.0)

    # Merge bin-level x back to the exploded table
    merge_cols = list(key_cols)
    out = out.merge(
        g[merge_cols + ["x_coord", "bin_width", "type", "is_highlight_bin"]],
        on=merge_cols,
        how="left",
        validate="many_to_one",
    )

    return out, g


def _compute_variable_x_on_bins(
    bins: pd.DataFrame,
    *,
    neutral_target_factor: float,
    backbone_factor: float,
    backbone_label: str = "backbone",
    genes_col: str = "genes",
) -> pd.DataFrame:
    """
    Compute bin_width/x_coord/type on a collapsed bin table where `genes` is a list per bin,
    and `is_highlight_bin` already exists.
    """
    out = bins.copy()

    def is_backbone(glist) -> bool:
        return isinstance(glist, list) and (backbone_label in glist)

    out["is_backbone_bin"] = out[genes_col].apply(is_backbone)
    out["type"] = np.where(out["is_backbone_bin"], "Backbone", "Target")

    width = np.full(len(out), float(neutral_target_factor), dtype=float)
    width[out["is_backbone_bin"].to_numpy()] = float(backbone_factor)
    width[out["is_highlight_bin"].to_numpy()] = 1.0

    out["bin_width"] = width
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

    Assumption: for a given bin, numeric values like log2/depth are identical across exploded rows.
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


def _add_cnvkit_segment_stub_legend(ax, segs_chr: pd.DataFrame | None) -> None:
    """
    Add a compact bottom-right legend with colored CNVkit segment stubs.
    Only shows entries that actually exist on this chromosome.
    """

    if segs_chr is None or segs_chr.empty or "cnvkit_cnv_call" not in segs_chr.columns:
        return

    calls = (
        segs_chr["cnvkit_cnv_call"].astype("string").fillna("").str.strip().str.upper()
    )

    handles = []

    if calls.eq("AMPLIFICATION").any():
        handles.append(Line2D([0], [0], color="red", lw=2.0, label="CNVKIT:Amp"))

    if calls.eq("DELETION").any():
        handles.append(Line2D([0], [0], color="royalblue", lw=2.0, label="CNVKIT:Del"))

    # Neutral/other segments
    if (~calls.isin(["AMPLIFICATION", "DELETION"]) | calls.eq("")).any():
        handles.append(Line2D([0], [0], color="black", lw=2.0, label="CNVKIT:Neutral"))

    if not handles:
        return

    legend = ax.legend(
        handles=handles,
        loc="lower right",
        fontsize=7,
        frameon=True,
        borderpad=0.3,
        labelspacing=0.3,
        handlelength=1.4,
        handletextpad=0.6,
    )

    # Ensure it stays even if you later add another legend
    ax.add_artist(legend)


def _annotate_highlight_bins(
    bins: pd.DataFrame,
    *,
    highlighted_genes: Iterable[str],
    genes_col: str = "genes",
    out_col: str = "is_highlight_bin",
) -> pd.DataFrame:
    """
    Add boolean out_col indicating whether any gene in bins[genes_col] is in highlighted_genes.
    """
    hi = set(map(str, highlighted_genes))

    def has_hi(glist) -> bool:
        if not isinstance(glist, list) or not glist:
            return False
        # intersection
        return any((g in hi) for g in glist)

    out = bins.copy()
    out[out_col] = out[genes_col].apply(has_hi)
    return out


def plot_chromosomes(
    cnr_df: pd.DataFrame,
    vcf_path: Path,
    gdf: pd.DataFrame,
    outdir: Path,
    case_id: str,
    pon_df: pd.DataFrame | None = None,
    gchunk: pd.DataFrame | None = None,
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
      - PON CNV chunks from gene_chunk_df (if provided)
      - gene highlighting + labels
      - BAF from VCF
    """
    outdir.mkdir(parents=True, exist_ok=True)

    # Rename columns for plotting:
    targets_col = "n.targets"

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

    # Compute gene flags & gene-level highlight decisions
    gdf = _compute_row_flags(gdf)
    gene_level = _compute_gene_level_highlights(
        gdf,
        targets_col=targets_col,
        highlight_only_cancer=highlight_only_cancer,
        min_gene_targets=MIN_GENE_TARGETS,
        min_gene_targets_cancer=MIN_GENE_TARGETS_CANCER,
    )

    stable_color = _stable_gene_color_fn()

    y_clip = float(y_abs_max)
    y_lim_chr = (-y_clip, y_clip)

    for chr_name in chr_order:
        # collapsed bin-level table for this chromosome
        sub_bins = merged[merged["chr"] == chr_name].copy()
        if sub_bins.empty:
            continue
        sub_bins = sub_bins.sort_values("start", kind="stable")

        g_chr = gdf[gdf["chr"] == chr_name].copy()
        g_chunks_chr = (
            gchunk[gchunk["chr"] == chr_name].copy()
            if (gchunk is not None and not gchunk.empty)
            else pd.DataFrame()
        )

        # --- Determine highlighted genes (unchanged logic) ---
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
        highlighted = highlighted[highlighted != "backbone"]
        gene_to_color = _make_gene_colors(highlighted)

        # Mark which collapsed bins contain any highlighted gene
        sub_bins = _annotate_highlight_bins(
            sub_bins,
            highlighted_genes=highlighted,
            genes_col="genes",
            out_col="is_highlight_bin",
        )

        # Compute variable x-coordinates on collapsed bins
        bins = _compute_variable_x_on_bins(
            sub_bins,
            neutral_target_factor=neutral_target_factor,
            backbone_factor=backbone_factor,
            backbone_label="backbone",
            genes_col="genes",
        )
        if bins.empty:
            continue

        # Smooth + clip on UNIQUE bins
        bins = bins.sort_values("x_coord", kind="stable").copy()
        bins["log2_smooth"] = bins["log2"].rolling(window=window, center=True).median()
        bins["log2_clipped"] = bins["log2"].clip(-y_clip, y_clip)
        bins["log2_smooth_clipped"] = bins["log2_smooth"].clip(-y_clip, y_clip)

        # Map function should use UNIQUE bins (stable mapping)
        pos_to_xcoord = _pos_to_xcoord_fn(bins)

        # VAF chr
        vaf_chr = vcf[vcf["CHROM"] == chr_name].sort_values("POS").copy()
        if not vaf_chr.empty:
            vaf_chr["x_coord"] = vaf_chr["POS"].apply(lambda p: pos_to_xcoord(int(p)))

        # Segments from g_chr
        segs_chr = _collect_segments_for_chr(g_chr, chr_name, y_clip=y_clip)
        if not segs_chr.empty:
            segs_chr = segs_chr.copy()
            segs_chr["x_start"] = segs_chr["cnvkit_seg_start"].apply(
                lambda p: pos_to_xcoord(int(p))
            )
            segs_chr["x_end"] = segs_chr["cnvkit_seg_end"].apply(
                lambda p: pos_to_xcoord(int(p))
            )

        # --- figure ---
        fig, (ax1, ax2) = plt.subplots(
            2, 1, figsize=(14, 6), sharex=True, gridspec_kw={"height_ratios": [2, 1]}
        )

        x = bins["x_coord"]

        # background should be BIN-LEVEL and exclude bins that have any highlighted gene
        bg_bins = bins[~bins["is_highlight_bin"]].copy()
        _draw_background_bins(ax1, bg_bins, y_col="log2_clipped")

        # highlighted bins: one point per BIN, gene-colored if unique highlighted gene else neutral
        if highlighted.size > 0:
            _draw_highlighted_bins(
                ax1,
                bins=bins,
                highlighted_genes=highlighted,
                gene_to_color=gene_to_color,
                genes_col="genes",
                y_col_bins="log2_clipped",
            )

        y_min, y_max = y_lim_chr

        if use_pon:
            _draw_pon_bars(ax1, bins, x, y_clip=y_clip)

        ax1.plot(
            x,
            bins["log2_smooth_clipped"],
            linewidth=1.5,
            alpha=0.9,
            color="tab:green",
            label=f"log2 (median {window} bins)",
        )

        pon_cnv_genes, _ = _draw_chunk_pon_segments(
            ax1,
            g_chunks_chr=g_chunks_chr,
            sub=bins,  # span_start/end should be based on bins
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
        ax1.set_title(title + "\n")
        ax1.legend(loc="upper right", fontsize=8)

        label_genes = sorted(set(map(str, highlighted)))
        label_genes = [g for g in label_genes if g != "backbone"]
        if label_genes:
            _draw_gene_labels(
                ax1,
                g_chr=g_chr,
                label_genes=label_genes,
                gene_to_color=gene_to_color,
                stable_color_fn=stable_color,
                pos_to_xcoord=pos_to_xcoord,
                y_max=y_max,
                base_label_offset=base_label_offset,
            )

        _draw_segments(ax1, segs_chr)
        _add_cnvkit_segment_stub_legend(ax1, segs_chr)
        _plot_vaf_panel(ax2, vaf_chr)
        ax2.set_xlabel(
            "Pseudo-position (highlighted genes expanded, other bins compressed)"
        )

        plt.tight_layout()
        out_png = outdir / f"cnv_chr{chr_name}_segments.png"
        plt.savefig(out_png, dpi=150)
        plt.close(fig)
