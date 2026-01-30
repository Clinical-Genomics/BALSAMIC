from __future__ import annotations

# Standard library
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

# Third-party
import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError
import fitz

# Local
from BALSAMIC.constants.analysis import Gender


FINAL_FLOAT_COLUMNS = [
    "seg_baf",
    "mean_log2",
    "min_log2",
    "max_log2",
    "mean_weight",
    "pon_gene_mean_log2",
    "pon_gene_mean_spread",
    "pon_gene_effect",
    "pon_gene_z",
    "depth_mean",
    "pon_mean_log2",
    "pon_mean_spread",
    "pon_chunk_effect",
    "pon_chunk_z",
]


# =============================================================================
# Generic helpers
# =============================================================================


def safe_read_csv(path: str | Path, **kwargs) -> pd.DataFrame:
    """Read CSV/TSV; return empty DataFrame if file is empty or missing."""
    try:
        return pd.read_csv(path, **kwargs)
    except (EmptyDataError, FileNotFoundError):
        return pd.DataFrame()


def _strip_chr_prefix(series: pd.Series) -> pd.Series:
    """Normalize chromosome values by stripping a leading 'chr' prefix."""
    return series.astype(str).str.replace("^chr", "", regex=True)


def _chrom_sort_key(chrom: str) -> tuple[int, int | str]:
    """Stable sort key: autosomes numeric first, then X/Y, then other contigs."""
    try:
        return (0, int(chrom))
    except ValueError:
        c = str(chrom)
        if c in ("X", "x"):
            return (1, 23)
        if c in ("Y", "y"):
            return (1, 24)
        return (2, c)


def _stable_sort_by_chr_interval(
    df: pd.DataFrame,
    chr_col: str,
    start_col: str,
    end_col: str,
    *,
    tmp_col: str = "chr_sort",
) -> pd.DataFrame:
    """Stable-sort by (chr, start, end) using `_chrom_sort_key`."""
    df = df.copy()
    df[tmp_col] = df[chr_col].map(_chrom_sort_key)
    df = df.sort_values(by=[tmp_col, start_col, end_col], kind="stable").drop(
        columns=[tmp_col]
    )
    return df


def _reorder_columns(df: pd.DataFrame, preferred: list[str]) -> pd.DataFrame:
    """Move preferred columns first (when present), preserving the rest."""
    preferred_present = [c for c in preferred if c in df.columns]
    return df[preferred_present + [c for c in df.columns if c not in preferred_present]]


def _flatten_agg_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Flatten MultiIndex columns produced by pandas .agg()."""
    df = df.copy()
    df.columns = [
        "_".join(col).strip("_") if isinstance(col, tuple) else col
        for col in df.columns
    ]
    return df


def _pon_abs_z(effect: float, spread: float, n_targets: float, *, min_n: int) -> float:
    """Compute abs(effect) / (spread / sqrt(n)) with defensive checks."""
    if pd.isna(effect) or pd.isna(spread) or spread <= 0:
        return np.nan
    if pd.isna(n_targets) or n_targets < min_n:
        return np.nan
    sigma_eff = spread / np.sqrt(float(n_targets))
    if sigma_eff <= 0:
        return np.nan
    return abs(effect) / sigma_eff


def _pon_direction(effect: float) -> str:
    """Map effect sign to 'gain'/'loss'/'neutral' (or '' when missing)."""
    if pd.isna(effect):
        return ""
    if effect > 0:
        return "gain"
    if effect < 0:
        return "loss"
    return "neutral"


def _pon_significance(z: float, *, noise_lt: float, borderline_lt: float) -> str:
    """Map z-score to 'noise'/'borderline'/'significant' (or '' when missing)."""
    if pd.isna(z):
        return ""
    if z < noise_lt:
        return "noise"
    if z < borderline_lt:
        return "borderline"
    return "significant"


def _pon_cnv_call_from_log2(
    *,
    is_significant: bool,
    mean_log2: float,
    gain_gt: float,
    loss_lt: float,
    non_significant_value: str,
    neutral_value: str,
) -> str:
    """Return AMPLIFICATION/DELETION/neutral/nonsignificant depending on thresholds."""
    if not is_significant or pd.isna(mean_log2):
        return non_significant_value
    if mean_log2 > gain_gt:
        return "AMPLIFICATION"
    if mean_log2 < loss_lt:
        return "DELETION"
    return neutral_value


# =============================================================================
# Shared IO / normalization (removes duplication)
# =============================================================================


@dataclass(frozen=True)
class CnrMeta:
    has_depth: bool
    has_weight: bool


def _detect_chr_col(df: pd.DataFrame, candidates: list[str]) -> str:
    for c in candidates:
        if c in df.columns:
            return c
    raise ValueError(f"Could not find chromosome column among: {candidates}")


def _explode_multigene_bins(
    cnr: pd.DataFrame, *, gene_col: str = "gene"
) -> pd.DataFrame:
    """
    Drop pseudo genes and explode multi-gene bins into one row per gene.symbol.
    """
    cnr = cnr.copy()
    cnr = cnr[cnr[gene_col].notna()]
    cnr = cnr[~cnr[gene_col].isin(["Antitarget", "-"])]

    cnr["gene.symbol"] = (
        cnr[gene_col]
        .astype(str)
        .str.split(",")
        .apply(lambda xs: [g.strip() for g in xs if g and g.strip()])
    )
    cnr = cnr.explode("gene.symbol")
    return cnr


def _load_cnr_bins(cnr_path: str | Path) -> tuple[pd.DataFrame, CnrMeta]:
    """
    Load CNVkit CNR, normalize chr, validate required columns, explode genes.

    Returns:
      bins_df with canonical 'chr' column and 'gene.symbol' column,
      plus metadata flags (depth/weight availability).
    """
    cnr = safe_read_csv(cnr_path, sep="\t").copy()
    if cnr.empty:
        return pd.DataFrame(), CnrMeta(False, False)

    chr_col = _detect_chr_col(cnr, ["chromosome", "chr", "CHROM"])
    cnr[chr_col] = _strip_chr_prefix(cnr[chr_col])

    required = {"start", "end", "log2", "gene"}
    missing = required - set(cnr.columns)
    if missing:
        raise ValueError(f"CNR missing required columns: {missing}")

    meta = CnrMeta(
        has_depth=("depth" in cnr.columns), has_weight=("weight" in cnr.columns)
    )

    cnr = _explode_multigene_bins(cnr, gene_col="gene")
    cnr = cnr.rename(columns={chr_col: "chr"})
    return cnr, meta


def _load_pon_bins(pon_path: str | Path) -> pd.DataFrame:
    """
    Load CNVkit PON .cnn, normalize chr, validate required columns,
    rename log2/spread -> pon_log2/pon_spread.
    """
    pon = safe_read_csv(pon_path, sep="\t").copy()
    if pon.empty:
        return pd.DataFrame()

    chr_col = _detect_chr_col(pon, ["chromosome", "chr", "CHROM"])
    pon[chr_col] = _strip_chr_prefix(pon[chr_col])

    required = {"start", "end", "log2", "spread"}
    missing = required - set(pon.columns)
    if missing:
        raise ValueError(
            "PON file must contain 'start','end','log2','spread' columns, "
            f"missing: {missing}"
        )

    pon = pon.rename(
        columns={chr_col: "chr", "log2": "pon_log2", "spread": "pon_spread"}
    )
    return pon[["chr", "start", "end", "pon_log2", "pon_spread"]]


def _left_merge_pon(cnr_bins: pd.DataFrame, pon_bins: pd.DataFrame) -> pd.DataFrame:
    """Left-merge PON columns onto bins by (chr,start,end)."""
    if pon_bins is None or pon_bins.empty:
        return cnr_bins
    return cnr_bins.merge(pon_bins, how="left", on=["chr", "start", "end"])


def _load_cns_segments(cns_path: str | Path) -> pd.DataFrame:
    """
    Load CNVkit CNS, normalize chr, keep available segment columns,
    and assign stable segment_id within each file.
    """
    cns = safe_read_csv(cns_path, sep="\t").copy()
    if cns.empty:
        return cns

    chr_col = _detect_chr_col(cns, ["chromosome", "chr", "CHROM"])
    cns[chr_col] = _strip_chr_prefix(cns[chr_col])

    if not {"start", "end"}.issubset(cns.columns):
        raise ValueError("CNS file must contain 'start' and 'end' columns.")

    base_cols = ["start", "end", "log2", "baf", "cn", "cn1", "cn2", "depth"]
    keep = [c for c in base_cols if c in cns.columns]
    cns = cns[[chr_col] + keep].rename(columns={chr_col: "chr"}).copy()
    cns = cns.sort_values(["chr", "start"], kind="stable").reset_index(drop=True)
    cns["segment_id"] = cns.index.astype(str)
    return cns


# =============================================================================
# Vectorized bin->segment assignment (replaces per-row scan)
# =============================================================================

_SEG_OUT_COLS = [
    "segment_id",
    "seg_start",
    "seg_end",
    "seg_log2",
    "seg_baf",
    "seg_cn",
    "seg_cn1",
    "seg_cn2",
]


def _assign_segments_by_center(
    bins_chr: pd.DataFrame, segs_chr: pd.DataFrame
) -> pd.DataFrame:
    """
    Assign a segment to each bin by bin-center, for a single chromosome.

    - Uses searchsorted on seg start positions.
    - Expected segments are non-overlapping and sorted by start.
    """
    out = bins_chr.copy()

    if segs_chr is None or segs_chr.empty or out.empty:
        out["segment_id"] = "no_segment"
        for col in _SEG_OUT_COLS[1:]:
            out[col] = np.nan
        return out

    segs = segs_chr.sort_values("start", kind="stable").reset_index(drop=True)

    centers = ((out["start"].to_numpy() + out["end"].to_numpy()) // 2).astype(int)

    starts = segs["start"].to_numpy()
    ends = segs["end"].to_numpy()

    # Find candidate segment index: rightmost segment with start <= center
    idx = np.searchsorted(starts, centers, side="right") - 1
    valid = (idx >= 0) & (centers < ends[np.clip(idx, 0, len(ends) - 1)])

    # Defaults
    out["segment_id"] = "no_segment"
    out["seg_start"] = np.nan
    out["seg_end"] = np.nan
    out["seg_log2"] = np.nan
    out["seg_baf"] = np.nan
    out["seg_cn"] = np.nan
    out["seg_cn1"] = np.nan
    out["seg_cn2"] = np.nan

    if not valid.any():
        return out

    vidx = idx[valid]
    out.loc[valid, "segment_id"] = segs.loc[vidx, "segment_id"].to_numpy()
    out.loc[valid, "seg_start"] = segs.loc[vidx, "start"].to_numpy()
    out.loc[valid, "seg_end"] = segs.loc[vidx, "end"].to_numpy()

    if "log2" in segs.columns:
        out.loc[valid, "seg_log2"] = segs.loc[vidx, "log2"].to_numpy()
    if "baf" in segs.columns:
        out.loc[valid, "seg_baf"] = segs.loc[vidx, "baf"].to_numpy()
    if "cn" in segs.columns:
        out.loc[valid, "seg_cn"] = segs.loc[vidx, "cn"].to_numpy()
    if "cn1" in segs.columns:
        out.loc[valid, "seg_cn1"] = segs.loc[vidx, "cn1"].to_numpy()
    if "cn2" in segs.columns:
        out.loc[valid, "seg_cn2"] = segs.loc[vidx, "cn2"].to_numpy()

    return out


def _assign_segments_all_chrom(bins: pd.DataFrame, segs: pd.DataFrame) -> pd.DataFrame:
    """
    Apply bin->segment assignment per chromosome; keeps original bins order stable.
    """
    if bins.empty:
        return bins

    if segs is None or segs.empty:
        out = bins.copy()
        out["segment_id"] = "no_segment"
        for col in _SEG_OUT_COLS[1:]:
            out[col] = np.nan
        return out

    segs_by_chr = {c: d for c, d in segs.groupby("chr", sort=False)}

    annotated = []
    for chrom, bins_chr in bins.groupby("chr", sort=False):
        annotated.append(
            _assign_segments_by_center(bins_chr, segs_by_chr.get(chrom, pd.DataFrame()))
        )
    return pd.concat(annotated, ignore_index=True)


# =============================================================================
# CNV classification
# =============================================================================


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

    chrom_str = re.sub(r"^chr", "", str(chrom), flags=re.IGNORECASE)

    if chrom_str.isdigit():
        expected = 2
    elif chrom_str.upper() == "X":
        expected = 2 if sex == Gender.FEMALE else 1
    elif chrom_str.upper() == "Y":
        expected = 1 if sex == Gender.MALE else 0
    else:
        expected = 2

    cn_int = int(round(float(cn)))
    exp_int = int(round(expected))

    if cn_int < exp_int:
        return "DELETION"
    if cn_int > exp_int:
        return "AMPLIFICATION"
    return "NEUTRAL"


# =============================================================================
# refGene / exon coverage (unchanged API; just reused)
# =============================================================================


def load_refgene_exons(
    refgene_path: str | Path,
    transcript_selection: str = "longest_tx",
) -> Dict[Tuple[str, str], dict]:
    cols = [
        "gene_symbol",
        "transcript",
        "chrom",
        "strand",
        "txStart",
        "txEnd",
        "cdsStart",
        "cdsEnd",
        "exonCount",
        "exonStarts",
        "exonEnds",
    ]
    rg = pd.read_csv(refgene_path, sep="\t", header=None, names=cols)

    rg["chrom"] = rg["chrom"].astype(str).str.replace("^chr", "", regex=True)
    rg["txStart"] = rg["txStart"].astype(int)
    rg["txEnd"] = rg["txEnd"].astype(int)

    exon_map: Dict[Tuple[str, str], dict] = {}

    for (chrom, gene), sub in rg.groupby(["chrom", "gene_symbol"], sort=False):
        if transcript_selection == "longest_tx":
            sub = sub.copy()
            sub["tx_len"] = sub["txEnd"] - sub["txStart"]
            best = sub.loc[sub["tx_len"].idxmax()]
        else:
            best = sub.iloc[0]

        starts = [int(x) for x in str(best["exonStarts"]).split(",") if x != ""]
        ends = [int(x) for x in str(best["exonEnds"]).split(",") if x != ""]
        exons: List[Tuple[int, int]] = list(zip(starts, ends))
        if not exons:
            continue

        exons.sort(key=lambda x: x[0])
        if best["strand"] == "-":
            exons = exons[::-1]

        exon_map[(chrom, gene)] = {
            "transcript": best["transcript"],
            "strand": best["strand"],
            "exons": exons,
        }

    return exon_map


# =============================================================================
# Cytoband (kept as-is except tiny reuse of _strip_chr_prefix)
# =============================================================================


def load_cytobands(path: str | Path) -> pd.DataFrame:
    """Load cytoband UCSC file; return normalized df with integer coords."""
    cols = ["chrom", "chromStart", "chromEnd", "name", "gieStain"]
    cyto = pd.read_csv(path, sep="\t", header=None, names=cols)
    cyto["chrom"] = _strip_chr_prefix(cyto["chrom"])
    cyto["start_int"] = cyto["chromStart"].astype(int)
    cyto["end_int"] = cyto["chromEnd"].astype(int)
    return cyto


def annotate_genes_with_cytoband(
    df_genes: pd.DataFrame, cyto: pd.DataFrame
) -> pd.DataFrame:
    """
    Add 'cytoband' based on (seg_start, seg_end), fallback to (start, end).
    """
    df = df_genes.copy()
    if "chr" not in df.columns:
        raise ValueError("df_genes must contain a 'chr' column.")
    df["chr"] = _strip_chr_prefix(df["chr"])
    df["cytoband"] = pd.Series(pd.NA, index=df.index, dtype="string")
    if "seg_start" not in df.columns:
        df["seg_start"] = np.nan
    if "seg_end" not in df.columns:
        df["seg_end"] = np.nan

    for chrom, genes_chr in df.groupby("chr", sort=False):
        cyto_chr = cyto[cyto["chrom"] == chrom]
        if cyto_chr.empty:
            continue

        cyto_chr = cyto_chr.sort_values("start_int", kind="stable")
        band_starts = cyto_chr["start_int"].to_numpy()
        band_ends = cyto_chr["end_int"].to_numpy()
        band_names = cyto_chr["name"].to_numpy()

        for idx, row in genes_chr.iterrows():
            if pd.notna(row["seg_start"]) and pd.notna(row["seg_end"]):
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
                label = f"{chrom}{names[0]}-{chrom}{names[-1]}"
            df.at[idx, "cytoband"] = label

    return df


def pdf_first_page_to_png(
    pdf_path: str | Path, png_path: str | Path, dpi: int = 300
) -> None:
    """
    Render the first page of a PDF to a PNG image using PyMuPDF.
    """
    doc = fitz.open(str(pdf_path))
    try:
        page = doc[0]
        page.get_pixmap(dpi=dpi).save(str(png_path))
    finally:
        doc.close()


def load_cancer_gene_set(
    path: str | Path, min_occurrence: int = 1, only_annotated: bool = True
) -> set[str]:
    """
    Load and filter a cancer gene list TSV into a set of gene symbols.

    Applies an occurrence threshold, optional OncoKB annotation filter,
    and optional restriction to classic cancer gene types (ONCOGENE / TSG)
    when the relevant columns exist.

    Parameters
    ----------
    path : str | Path
        Path to cancer gene list TSV file.
    min_occurrence : int, optional
        Minimum occurrence count required to include a gene (default: 1).
    only_annotated : bool, optional
        If True, require OncoKB annotation = YES when the column is present
        (default: True).

    Returns
    -------
    set[str]
        Gene symbols passing filters.
    """
    df = pd.read_csv(path, sep="\t", dtype=str)
    df.columns = df.columns.str.strip()

    symbol_col = "Hugo Symbol"
    occ_col = "# of occurrence within resources (Column J-P)"
    onco_col = "OncoKB Annotated"
    type_col = "Gene Type"

    if symbol_col not in df.columns or occ_col not in df.columns:
        missing = [c for c in (symbol_col, occ_col) if c not in df.columns]
        raise ValueError(f"Cancer gene TSV missing required columns: {missing}")

    df[occ_col] = pd.to_numeric(df[occ_col], errors="coerce").fillna(0).astype(int)

    mask = df[occ_col] >= int(min_occurrence)

    if only_annotated and onco_col in df.columns:
        mask &= df[onco_col].astype(str).str.upper().eq("YES")

    if type_col in df.columns:
        mask &= df[type_col].isin(["ONCOGENE", "TSG"])

    selected = df.loc[mask, symbol_col].astype(str).str.strip()
    selected = selected[selected != ""]
    return set(selected)


def compute_dlr_spread_from_cnr(
    cnr_path: str | Path,
    weight_thresh: float = 0.1,
    targets_only: bool = True,
    exclude_sex_chromosomes: bool = True,
) -> float:
    """
    Compute a DLR-like spread metric from a CNVkit .cnr file.

    - Sort bins by chr, start.
    - Optionally keep only target bins, with weight > threshold, autosomes only.
    - Compute differences between adjacent log2 values.
    - Return np.nanstd of these differences.

    Returns np.nan if fewer than 3 usable bins.
    """
    cnr = pd.read_csv(cnr_path, sep="\t")
    if cnr.empty or "log2" not in cnr.columns:
        return float("nan")

    chr_col = "chromosome" if "chromosome" in cnr.columns else "chr"
    cnr[chr_col] = _strip_chr_prefix(cnr[chr_col])

    # mark Target / Antitarget
    if "gene" in cnr.columns:
        cnr["type"] = np.where(cnr["gene"] == "Antitarget", "Antitarget", "Target")
    else:
        cnr["type"] = "Target"

    if targets_only:
        cnr = cnr[cnr["type"] == "Target"]

    if "weight" in cnr.columns:
        cnr = cnr[cnr["weight"] > weight_thresh]

    if exclude_sex_chromosomes:
        cnr = cnr[~cnr[chr_col].isin(["X", "x", "Y", "y"])]

    cnr = cnr.dropna(subset=["log2"])
    if cnr.empty:
        return float("nan")

    cnr["chr_sort"] = cnr[chr_col].map(_chrom_sort_key)
    cnr = cnr.sort_values(["chr_sort", "start"]).reset_index(drop=True)
    log2_vals = cnr["log2"].to_numpy()

    if log2_vals.size < 3:
        return float("nan")

    diffs = np.diff(log2_vals)
    dlr = float(np.nanstd(diffs))
    return dlr


def compute_pon_spread_summaries(
    cnn_path: str | Path,
    targets_only: bool = True,
    exclude_sex_chromosomes: bool = True,
) -> dict[str, float]:
    """
    Compute simple spread summaries from a CNVkit PON .cnn file.

    Returns a dict with:
      - pon_spread_median_target
      - pon_spread_q90_target

    (np.nan if not available)
    """
    cnn = pd.read_csv(cnn_path, sep="\t")
    if cnn.empty or "spread" not in cnn.columns:
        return {
            "pon_spread_median_target": float("nan"),
            "pon_spread_q90_target": float("nan"),
        }

    chr_col = "chromosome" if "chromosome" in cnn.columns else "chr"
    cnn[chr_col] = cnn[chr_col].astype(str).str.replace("^chr", "", regex=True)

    # mark Target / Antitarget
    if "gene" in cnn.columns:
        cnn["type"] = np.where(cnn["gene"] == "Antitarget", "Antitarget", "Target")
    else:
        cnn["type"] = "Target"

    if targets_only:
        cnn = cnn[cnn["type"] == "Target"]

    if exclude_sex_chromosomes:
        cnn = cnn[~cnn[chr_col].isin(["X", "x", "Y", "y"])]

    cnn = cnn.dropna(subset=["spread"])
    if cnn.empty:
        return {
            "pon_spread_median_target": float("nan"),
            "pon_spread_q90_target": float("nan"),
        }

    return {
        "pon_spread_median_target": float(cnn["spread"].median()),
        "pon_spread_q90_target": float(cnn["spread"].quantile(0.90)),
    }


def compute_summary_metrics(
    cnr_path: str | Path,
    cnn_path: str | Path | None,
    gene_seg_df: pd.DataFrame | None,
) -> pd.DataFrame:
    """
    Build a 1-row DataFrame with QC / burden summaries.

    Includes:
      - Log2-spread (DLR-like) from the CNVkit .cnr file
      - PON spread summaries from the CNVkit .cnn file (targets only)
      - Gene-level CNV/LOH counts from the gene-segment table (if provided)

    Parameters
    ----------
    cnr_path : str | Path
        Path to CNVkit .cnr file.
    cnn_path : str | Path | None
        Path to CNVkit PON .cnn file (optional).
    gene_seg_df : pd.DataFrame | None
        Gene-segment table (optional).

    Returns
    -------
    pd.DataFrame
        Single-row DataFrame of metrics.
    """
    metrics: dict[str, float | int] = {}

    # 1) DLR-like spread
    metrics["Log2-spread"] = compute_dlr_spread_from_cnr(cnr_path)

    # 2) PON spread summaries
    if cnn_path is not None and Path(cnn_path).is_file():
        pon_stats = compute_pon_spread_summaries(cnn_path)
        metrics.update(pon_stats)
    else:
        metrics["pon_spread_median_target"] = float("nan")
        metrics["pon_spread_q90_target"] = float("nan")

    # 3) Gene-level CNV / LOH summaries
    gene_stats = compute_gene_cnv_summaries(gene_seg_df)
    metrics.update(gene_stats)

    return pd.DataFrame([metrics])


# =============================================================================
# Gene table helpers (shared between gene-seg and gene-chunk)
# =============================================================================


def _compute_gene_level_pon_from_bins(bins: pd.DataFrame) -> pd.DataFrame | None:
    """
    Gene-level PON stats computed from per-bin data (independent of segments/chunks).
    """
    if not {"pon_log2", "pon_spread"}.issubset(bins.columns):
        return None

    has_pon = bins["pon_spread"].notna()
    if not has_pon.any():
        return None

    gene_agg = (
        bins.loc[has_pon]
        .groupby(["chr", "gene.symbol"], as_index=False)
        .agg(
            gene_n_targets=("log2", "count"),
            gene_mean_log2=("log2", "mean"),
            pon_gene_mean_log2=("pon_log2", "mean"),
            pon_gene_mean_spread=("pon_spread", "mean"),
        )
    )

    MIN_N_FOR_PON = 2
    gene_agg["pon_gene_effect"] = (
        gene_agg["gene_mean_log2"] - gene_agg["pon_gene_mean_log2"]
    )

    gene_agg["pon_gene_z"] = gene_agg.apply(
        lambda r: _pon_abs_z(
            r["pon_gene_effect"],
            r["pon_gene_mean_spread"],
            r["gene_n_targets"],
            min_n=MIN_N_FOR_PON,
        ),
        axis=1,
    )

    gene_agg["pon_gene_direction"] = gene_agg["pon_gene_effect"].apply(_pon_direction)
    gene_agg["pon_gene_call"] = gene_agg["pon_gene_z"].apply(
        lambda z: _pon_significance(z, noise_lt=1.5, borderline_lt=3.0)
    )

    gene_agg["pon_gene_cnv_call"] = gene_agg.apply(
        lambda r: _pon_cnv_call_from_log2(
            is_significant=str(r.get("pon_gene_call", "")).strip().lower()
            == "significant",
            mean_log2=r.get("gene_mean_log2", np.nan),
            gain_gt=0.08,
            loss_lt=-0.08,
            non_significant_value="",
            neutral_value="",
        ),
        axis=1,
    )

    return gene_agg[
        [
            "chr",
            "gene.symbol",
            "pon_gene_mean_log2",
            "pon_gene_mean_spread",
            "pon_gene_effect",
            "pon_gene_z",
            "pon_gene_direction",
            "pon_gene_call",
            "pon_gene_cnv_call",
        ]
    ]


def _add_cnv_calls_from_total_cn(genes_df: pd.DataFrame, sex: Gender) -> pd.DataFrame:
    """
    Add cnvkit_cnv_call and purecn_cnv_call with the same "AMP_EPS" sanity check logic.
    """
    AMP_EPS = 0.05
    out = genes_df.copy()

    def _call_with_sanity(total_cn_col: str, sanity_log2_col: str) -> pd.Series:
        if not {"chr", total_cn_col}.issubset(out.columns):
            return pd.Series(index=out.index, dtype="object")

        def _row_call(r: pd.Series) -> str:
            base = classify_cnv_from_total_cn_sex_aware(r[total_cn_col], r["chr"], sex)
            base_u = str(base).strip().upper()
            if base_u not in ("AMPLIFICATION", "DELETION"):
                return base
            v = r.get(sanity_log2_col, np.nan)
            if pd.notna(v) and abs(float(v)) < AMP_EPS:
                return "NEUTRAL"
            return base

        return out.apply(_row_call, axis=1)

    if {"seg_cn", "chr"}.issubset(out.columns):
        out["cnvkit_cnv_call"] = _call_with_sanity("seg_cn", "seg_log2")

    if {"loh_C", "chr"}.issubset(out.columns):
        out["purecn_cnv_call"] = _call_with_sanity("loh_C", "loh_seg_mean")

    return out


def _finalize_gene_table_order(genes_df: pd.DataFrame) -> pd.DataFrame:
    """
    Shared final ordering/sorting.
    """
    out = genes_df.copy()
    if "n_targets" in out.columns and "n.targets" not in out.columns:
        out = out.rename(columns={"n_targets": "n.targets"})

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
        "cnvkit_cnv_call",
        "purecn_cnv_call",
        "pon_gene_mean_log2",
        "pon_gene_mean_spread",
        "pon_gene_effect",
        "pon_gene_z",
        "pon_gene_direction",
        "pon_gene_call",
        "pon_gene_cnv_call",
        "pon_mean_log2",
        "pon_mean_spread",
        "pon_chunk_effect",
        "pon_chunk_z",
    ]
    out = _reorder_columns(out, cols_order)
    return _stable_sort_by_chr_interval(out, "chr", "region_start", "region_end")


# =============================================================================
# build_gene_segment_table (refactored)
# =============================================================================


def build_gene_segment_table(
    cnr_path: str | Path,
    cns_path: str | Path,
    cancer_genes: Optional[Set[str]] = None,
    refgene_path: str | Path | None = None,
    transcript_selection: str = "longest_tx",
    loh_path: str | Path | None = None,
    cytoband_path: str | Path | None = None,
    sex: Gender = None,
    pon_path: str | Path | None = None,
) -> pd.DataFrame:
    """
    Build a per-gene × segment table from CNVkit CNR + CNS
    (+ optional refGene, LOHgenes, cytobands, PON).
    """

    # 1) CNR bins (+ explode genes)
    bins, meta = _load_cnr_bins(cnr_path)
    if bins.empty:
        return pd.DataFrame()

    # 2) optional PON -> bins
    if pon_path is not None and Path(pon_path).is_file():
        pon = _load_pon_bins(pon_path)
        bins = _left_merge_pon(bins, pon)

    # 3) CNS segments (optional)
    segs = _load_cns_segments(cns_path)

    # If CNS missing -> one row per gene over all bins
    if segs.empty:
        agg = {
            "region_start": ("start", "min"),
            "region_end": ("end", "max"),
            "n_targets": ("log2", "count"),
            "mean_log2": ("log2", "mean"),
            "min_log2": ("log2", "min"),
            "max_log2": ("log2", "max"),
        }
        if meta.has_weight and "weight" in bins.columns:
            agg["mean_weight"] = ("weight", "mean")

        genes_df = bins.groupby(["chr", "gene.symbol"], as_index=False).agg(**agg)

        for col in [
            "seg_start",
            "seg_end",
            "seg_log2",
            "seg_baf",
            "seg_cn",
            "seg_cn1",
            "seg_cn2",
        ]:
            genes_df[col] = np.nan

        if cancer_genes is not None:
            genes_df["is_cancer_gene"] = genes_df["gene.symbol"].isin(cancer_genes)

        # cytoband (optional)
        if cytoband_path is not None and Path(cytoband_path).is_file():
            cyto = load_cytobands(cytoband_path)
            genes_df = annotate_genes_with_cytoband(genes_df, cyto)

        return _finalize_gene_table_order(genes_df)

    # 4) assign segment to bins (vectorized)
    bins = _assign_segments_all_chrom(bins, segs)
    bins = bins.sort_values(["chr", "gene.symbol", "start"], kind="stable").reset_index(
        drop=True
    )

    # 5) gene-level PON stats (independent of segments)
    gene_pon = _compute_gene_level_pon_from_bins(bins)

    # 6) collapse bins -> (chr, gene.symbol, segment)
    agg_dict: dict[str, list[str]] = {
        "start": ["min"],
        "end": ["max"],
        "log2": ["count", "mean", "min", "max"],
    }
    if meta.has_depth and "depth" in bins.columns:
        agg_dict["depth"] = ["mean"]
    if meta.has_weight and "weight" in bins.columns:
        agg_dict["weight"] = ["mean"]

    # segment fields as first()
    for col in [
        "seg_start",
        "seg_end",
        "seg_log2",
        "seg_baf",
        "seg_cn",
        "seg_cn1",
        "seg_cn2",
    ]:
        if col in bins.columns:
            agg_dict[col] = ["first"]

    genes_df = bins.groupby(["chr", "gene.symbol", "segment_id"], as_index=False).agg(
        agg_dict
    )
    genes_df = _flatten_agg_columns(genes_df).rename(
        columns={
            "start_min": "region_start",
            "end_max": "region_end",
            "log2_count": "n_targets",
            "log2_mean": "mean_log2",
            "log2_min": "min_log2",
            "log2_max": "max_log2",
            "weight_mean": "mean_weight",
            "seg_start_first": "seg_start",
            "seg_end_first": "seg_end",
            "seg_log2_first": "seg_log2",
            "seg_baf_first": "seg_baf",
            "seg_cn_first": "seg_cn",
            "seg_cn1_first": "seg_cn1",
            "seg_cn2_first": "seg_cn2",
            "depth_mean": "depth_mean",
        }
    )
    genes_df = genes_df.drop(columns=["segment_id"], errors="ignore")

    if cancer_genes is not None:
        genes_df["is_cancer_gene"] = genes_df["gene.symbol"].isin(cancer_genes)

    # 7) LOHgenes attach (kept identical behavior, just clearer)
    if loh_path is not None and Path(loh_path).is_file():
        loh = safe_read_csv(loh_path, sep=",").copy()
        if not loh.empty:
            loh_chr_col = "chr" if "chr" in loh.columns else "chromosome"
            loh[loh_chr_col] = _strip_chr_prefix(loh[loh_chr_col])
            loh = loh.rename(columns={loh_chr_col: "chr"})

            if "Sampleid" in loh.columns and loh["Sampleid"].nunique() > 1:
                raise ValueError(
                    "LOHgenes file contains multiple Sampleid values; "
                    "subset to a single sample before calling build_gene_segment_table."
                )

            required_loh = {
                "gene.symbol",
                "chr",
                "C",
                "M",
                "M.flagged",
                "loh",
                "seg.mean",
                "num.snps",
            }
            missing = required_loh - set(loh.columns)
            if missing:
                raise ValueError(f"LOHgenes file missing columns: {missing}")

            loh_small = loh[
                [
                    "chr",
                    "gene.symbol",
                    "C",
                    "M",
                    "M.flagged",
                    "loh",
                    "seg.mean",
                    "num.snps",
                ]
            ].rename(
                columns={
                    "C": "loh_C",
                    "M": "loh_M",
                    "M.flagged": "loh_M_flagged",
                    "loh": "loh_flag",
                    "seg.mean": "loh_seg_mean",
                    "num.snps": "loh_num_snps",
                }
            )
            genes_df = genes_df.merge(loh_small, how="left", on=["chr", "gene.symbol"])

    # 8) cytoband
    if cytoband_path is not None and Path(cytoband_path).is_file():
        cyto = load_cytobands(cytoband_path)
        genes_df = annotate_genes_with_cytoband(genes_df, cyto)

    # 9) CNV calls
    genes_df = _add_cnv_calls_from_total_cn(genes_df, sex)

    # 10) exon annotations (kept behavior: writes exons_hit)
    if refgene_path is not None and Path(refgene_path).is_file():
        exon_map = load_refgene_exons(
            refgene_path, transcript_selection=transcript_selection
        )

        def _segment_is_cnv(row: pd.Series) -> bool:
            seg_cn = row.get("seg_cn", np.nan)
            if pd.notna(seg_cn):
                call = classify_cnv_from_total_cn_sex_aware(seg_cn, row.get("chr"), sex)
                if str(call).upper() in ("AMPLIFICATION", "DELETION"):
                    return True
            loh_C = row.get("loh_C", np.nan)
            if pd.notna(loh_C):
                call = classify_cnv_from_total_cn_sex_aware(loh_C, row.get("chr"), sex)
                if str(call).upper() in ("AMPLIFICATION", "DELETION"):
                    return True
            return False

        def _exons_hit(row: pd.Series) -> str:
            key = (str(row["chr"]), str(row["gene.symbol"]))
            info = exon_map.get(key)
            if info is None:
                return ""

            use_segment_span = _segment_is_cnv(row)

            seg_start = row.get("seg_start")
            seg_end = row.get("seg_end")
            chunk_start = row.get("region_start")
            chunk_end = row.get("region_end")

            if use_segment_span and pd.notna(seg_start) and pd.notna(seg_end):
                region_start = seg_start
                region_end = seg_end
            else:
                region_start = chunk_start
                region_end = chunk_end

            if pd.isna(region_start) or pd.isna(region_end):
                return ""

            exons = info["exons"]
            gene_start = min(s for s, _e in exons)
            gene_end = max(e for _s, e in exons)

            if region_start <= gene_start and region_end >= gene_end:
                return "whole_gene"

            hit = [
                i
                for i, (s, e) in enumerate(exons, start=1)
                if e > region_start and s < region_end
            ]
            if not hit:
                return ""

            # compress consecutive indices
            ranges = []
            a = b = hit[0]
            for i in hit[1:]:
                if i == b + 1:
                    b = i
                else:
                    ranges.append((a, b))
                    a = b = i
            ranges.append((a, b))

            parts = [str(x) if x == y else f"{x}-{y}" for x, y in ranges]
            return ",".join(parts)

        genes_df["exons_hit"] = genes_df.apply(_exons_hit, axis=1)

    # 11) attach gene-level PON stats
    if gene_pon is not None and not gene_pon.empty:
        genes_df = genes_df.merge(gene_pon, how="left", on=["chr", "gene.symbol"])

    return _finalize_gene_table_order(genes_df)


# =============================================================================
# build_gene_chunk_table (refactored)
# =============================================================================


def build_gene_chunk_table(
    cnr_path: str | Path,
    gene_seg_df: pd.DataFrame,
    pon_path: str | Path | None = None,
) -> pd.DataFrame:
    """
    Build per-gene, per-chunk table using CNR + PON (required),
    then annotate chunks using existing gene-segment table by overlap.
    """

    # 0) PON required
    if pon_path is None or not Path(pon_path).is_file():
        return pd.DataFrame()

    # 1) CNR bins (exploded)
    bins, meta = _load_cnr_bins(cnr_path)
    if bins.empty:
        return pd.DataFrame()

    # 2) PON bins required + merge
    pon = _load_pon_bins(pon_path)
    if pon.empty:
        return pd.DataFrame()

    bins = _left_merge_pon(bins, pon)
    bins = bins.sort_values(["chr", "gene.symbol", "start"], kind="stable").reset_index(
        drop=True
    )

    # If all pon_spread missing -> unusable
    if "pon_spread" not in bins.columns or bins["pon_spread"].isna().all():
        return pd.DataFrame()

    # -------------------- 3) Within-gene chunking (same logic, less nesting) -------------------- #
    MIN_GENE_TARGETS = 8
    MIN_RUN_BINS = 4
    Z_BIN_THRESH = 1.5
    Z_RUN_THRESH = 3.0
    SMOOTH_WINDOW = 3

    bins["chunk_id"] = "no_chunk"
    next_id = 0

    def _run_z_for_indices(ix: list[int]) -> float:
        sub = bins.loc[ix]
        eff = sub["log2"].to_numpy() - sub["pon_log2"].fillna(0.0).to_numpy()
        if eff.size == 0:
            return np.nan
        mean_eff = float(np.nanmean(eff))
        if not np.isfinite(mean_eff):
            return np.nan

        sigma = sub["pon_spread"].fillna(0.0).to_numpy()
        sigma_safe = np.where(sigma <= 0, 1e-3, sigma)
        sigma_mean = float(np.nanmean(sigma_safe)) if sigma_safe.size else 1e-3

        n = len(eff)
        sigma_eff = sigma_mean / np.sqrt(float(n))
        if sigma_eff <= 0:
            return np.nan
        return abs(mean_eff) / sigma_eff

    def _assign_gene_chunks(df_gene: pd.DataFrame) -> None:
        nonlocal next_id

        if df_gene.shape[0] < MIN_GENE_TARGETS:
            return
        if df_gene["pon_spread"].isna().all():
            return

        idxs = df_gene.index.to_numpy()
        eff = df_gene["log2"].to_numpy() - df_gene["pon_log2"].fillna(0.0).to_numpy()
        sigma = df_gene["pon_spread"].fillna(0.0).to_numpy()
        sigma_safe = np.where(sigma <= 0, 1e-3, sigma)
        z_raw = eff / sigma_safe

        z_smooth = (
            pd.Series(z_raw, index=df_gene.index)
            .rolling(window=SMOOTH_WINDOW, center=True, min_periods=1)
            .median()
            .to_numpy()
        )

        run: list[int] = []
        run_sign: int | None = None

        def _finalize(ix: list[int]) -> None:
            nonlocal next_id
            if len(ix) < MIN_RUN_BINS:
                return
            rz = _run_z_for_indices(ix)
            if not np.isfinite(rz) or rz < Z_RUN_THRESH:
                return
            label = f"genechunk_{next_id}"
            next_id += 1
            bins.loc[ix, "chunk_id"] = label

        for i, z in zip(idxs, z_smooth):
            if not np.isfinite(z) or abs(z) < Z_BIN_THRESH:
                if run:
                    _finalize(run)
                    run = []
                    run_sign = None
                continue

            sgn = 1 if z > 0 else -1
            if run_sign is None or sgn == run_sign:
                run.append(i)
                run_sign = sgn
            else:
                _finalize(run)
                run = [i]
                run_sign = sgn

        if run:
            _finalize(run)

    for (_ch, _g), df_gene in bins.groupby(["chr", "gene.symbol"], sort=False):
        _assign_gene_chunks(df_gene)

    # -------------------- 3b) Bridge cleanup (same logic, extracted primitives) -------------------- #
    MAX_BRIDGE_BINS = 4
    BRIDGE_DELTA = 0.12
    SMALL_SEG_N = 3
    MERGE_DELTA = 0.10

    def _gene_runs(df_gene: pd.DataFrame) -> list[dict]:
        eff = df_gene["log2"].to_numpy() - df_gene["pon_log2"].fillna(0.0).to_numpy()
        labels = df_gene["chunk_id"].tolist()

        runs: list[dict] = []
        pos = [0]
        cur = labels[0]
        for i in range(1, len(labels)):
            if labels[i] == cur:
                pos.append(i)
            else:
                ix = df_gene.index[pos].tolist()
                mean_eff = float(np.nanmean(eff[pos])) if pos else np.nan
                runs.append({"positions": pos, "indices": ix, "mean_eff": mean_eff})
                cur = labels[i]
                pos = [i]
        if pos:
            ix = df_gene.index[pos].tolist()
            mean_eff = float(np.nanmean(eff[pos])) if pos else np.nan
            runs.append({"positions": pos, "indices": ix, "mean_eff": mean_eff})
        return runs, eff

    def _cleanup_gene_chunks(df_gene: pd.DataFrame) -> None:
        if df_gene.empty or df_gene.shape[0] == 1:
            return

        runs, eff_all = _gene_runs(df_gene)
        if len(runs) <= 1:
            return

        new_runs: list[dict] = []
        i = 0
        while i < len(runs):
            if i <= len(runs) - 3:
                A, B, C = runs[i], runs[i + 1], runs[i + 2]
                if (
                    len(B["indices"]) <= MAX_BRIDGE_BINS
                    and np.isfinite(A["mean_eff"])
                    and np.isfinite(C["mean_eff"])
                    and abs(A["mean_eff"] - C["mean_eff"]) <= BRIDGE_DELTA
                ):
                    merged_pos = A["positions"] + B["positions"] + C["positions"]
                    merged_ix = A["indices"] + B["indices"] + C["indices"]
                    merged_mean = (
                        float(np.nanmean(eff_all[merged_pos])) if merged_pos else np.nan
                    )
                    new_runs.append(
                        {
                            "positions": merged_pos,
                            "indices": merged_ix,
                            "mean_eff": merged_mean,
                        }
                    )
                    i += 3
                    continue

            r = runs[i]
            if new_runs:
                last = new_runs[-1]
                if (
                    np.isfinite(r["mean_eff"])
                    and np.isfinite(last["mean_eff"])
                    and abs(r["mean_eff"] - last["mean_eff"]) <= MERGE_DELTA
                    and (
                        len(r["indices"]) <= SMALL_SEG_N
                        or len(last["indices"]) <= SMALL_SEG_N
                    )
                ):
                    merged_pos = last["positions"] + r["positions"]
                    merged_ix = last["indices"] + r["indices"]
                    merged_mean = (
                        float(np.nanmean(eff_all[merged_pos])) if merged_pos else np.nan
                    )
                    new_runs[-1] = {
                        "positions": merged_pos,
                        "indices": merged_ix,
                        "mean_eff": merged_mean,
                    }
                else:
                    new_runs.append(r)
            else:
                new_runs.append(r)

            i += 1

        for run_idx, run in enumerate(new_runs):
            bins.loc[run["indices"], "chunk_id"] = f"genechunk_clean_{run_idx}"

    for (_ch, _g), df_gene in bins.groupby(["chr", "gene.symbol"], sort=False):
        _cleanup_gene_chunks(df_gene)

    bins = bins.sort_values(["chr", "gene.symbol", "start"], kind="stable").reset_index(
        drop=True
    )

    # -------------------- 4) Collapse to gene × chunk -------------------- #
    agg_dict: dict[str, list[str]] = {
        "start": ["min"],
        "end": ["max"],
        "log2": ["count", "mean", "min", "max"],
        "pon_log2": ["mean"],
        "pon_spread": ["mean"],
    }
    if meta.has_depth and "depth" in bins.columns:
        agg_dict["depth"] = ["mean"]
    if meta.has_weight and "weight" in bins.columns:
        agg_dict["weight"] = ["mean"]

    chunks_df = bins.groupby(["chr", "gene.symbol", "chunk_id"], as_index=False).agg(
        agg_dict
    )
    chunks_df = _flatten_agg_columns(chunks_df).rename(
        columns={
            "start_min": "region_start",
            "end_max": "region_end",
            "log2_count": "n_targets",
            "log2_mean": "mean_log2",
            "log2_min": "min_log2",
            "log2_max": "max_log2",
            "weight_mean": "mean_weight",
            "pon_log2_mean": "pon_mean_log2",
            "pon_spread_mean": "pon_mean_spread",
        }
    )
    chunks_df["n.targets"] = chunks_df["n_targets"]

    # -------------------- 5) PON deviation per chunk -------------------- #
    min_n_for_pon = 2
    chunks_df["pon_chunk_effect"] = chunks_df["mean_log2"] - chunks_df["pon_mean_log2"]
    chunks_df["pon_chunk_z"] = chunks_df.apply(
        lambda r: _pon_abs_z(
            r["pon_chunk_effect"],
            r["pon_mean_spread"],
            r.get("n_targets", r.get("n.targets", np.nan)),
            min_n=min_n_for_pon,
        ),
        axis=1,
    )
    chunks_df["pon_chunk_direction"] = chunks_df["pon_chunk_effect"].apply(
        _pon_direction
    )
    chunks_df["pon_chunk_significance"] = chunks_df["pon_chunk_z"].apply(
        lambda z: _pon_significance(z, noise_lt=2.0, borderline_lt=5.0)
    )
    chunks_df["pon_chunk_call"] = chunks_df.apply(
        lambda r: _pon_cnv_call_from_log2(
            is_significant=str(r.get("pon_chunk_significance", "")).strip().lower()
            == "significant",
            mean_log2=r.get("mean_log2", np.nan),
            gain_gt=0.07,
            loss_lt=-0.07,
            non_significant_value="NEUTRAL",
            neutral_value="NEUTRAL",
        ),
        axis=1,
    )

    # -------------------- 6) Annotate from gene_seg_df by overlap -------------------- #
    required = {"chr", "gene.symbol", "region_start", "region_end"}
    if not required.issubset(gene_seg_df.columns):
        raise ValueError(
            "gene_seg_df must contain 'chr', 'gene.symbol', 'region_start', 'region_end'."
        )

    gseg = gene_seg_df.copy()
    gseg["chr"] = _strip_chr_prefix(gseg["chr"])

    seg_map: dict[tuple[str, str], pd.DataFrame] = {
        (str(ch), str(g)): df.reset_index(drop=True)
        for (ch, g), df in gseg.groupby(["chr", "gene.symbol"], sort=False)
    }

    protected_cols = {
        "region_start",
        "region_end",
        "n_targets",
        "n.targets",
        "mean_log2",
        "min_log2",
        "max_log2",
        "mean_weight",
        "pon_mean_log2",
        "pon_mean_spread",
        "pon_chunk_effect",
        "pon_chunk_z",
        "pon_chunk_direction",
        "pon_chunk_significance",
        "pon_chunk_call",
    }

    annot_cols = [
        c
        for c in gseg.columns
        if c
        not in ({"chr", "gene.symbol", "region_start", "region_end"} | protected_cols)
    ]

    def _annotate_one(row: pd.Series) -> pd.Series:
        segs = seg_map.get((str(row["chr"]), str(row["gene.symbol"])))
        if segs is None or segs.empty:
            return row

        c0, c1 = row["region_start"], row["region_end"]
        if pd.isna(c0) or pd.isna(c1):
            return row

        mask = (segs["region_end"] > c0) & (segs["region_start"] < c1)
        if not mask.any():
            return row

        ssub = segs.loc[mask]
        ov0 = np.maximum(ssub["region_start"].to_numpy(), c0)
        ov1 = np.minimum(ssub["region_end"].to_numpy(), c1)
        best = int((ov1 - ov0).argmax())
        picked = ssub.iloc[best]

        for c in annot_cols:
            row[c] = picked[c]
        return row

    chunks_df = chunks_df.apply(_annotate_one, axis=1)

    # Drop any gene-level PON columns if they came in from gene_seg_df
    chunks_df = chunks_df.drop(
        columns=[
            "pon_gene_mean_log2",
            "pon_gene_mean_spread",
            "pon_gene_effect",
            "pon_gene_z",
            "pon_gene_direction",
            "pon_gene_call",
            "pon_gene_cnv_call",
            "pon_cnv_call",
        ],
        errors="ignore",
    )

    # -------------------- 7) Final sort/order -------------------- #
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
        "pon_chunk_effect",
        "pon_chunk_z",
        "pon_chunk_direction",
        "pon_chunk_significance",
        "pon_chunk_call",
        "cnvkit_cnv_call",
        "purecn_cnv_call",
    ]
    chunks_df = _reorder_columns(chunks_df, cols_order)
    chunks_df = _stable_sort_by_chr_interval(
        chunks_df, "chr", "region_start", "region_end"
    )
    chunks_df = chunks_df.drop(columns=["chunk_id", "n_targets"], errors="ignore")

    return chunks_df
