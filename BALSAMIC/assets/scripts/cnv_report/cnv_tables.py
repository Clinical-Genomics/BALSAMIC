from __future__ import annotations

# Standard library
from typing import Mapping

# Third-party
import numpy as np
import pandas as pd

# Local
from BALSAMIC.constants.analysis import Gender
from cnv_constants import (
    TableSpec,
    GENE_TABLE_SPEC,
    SEGMENT_TABLE_SPEC,
    CHR,
    GENE,
)
from cnv_report_utils import chrom_sort_key
from cnv_create_generegions import create_generegions


def finalize_table(df: pd.DataFrame, spec: TableSpec) -> pd.DataFrame:
    """
    Finalize a table according to a TableSpec.

    Applies:
      - preferred column ordering
      - float rounding
      - stable genomic sorting
    """
    out = df.copy()

    # Put preferred columns first, keep any extra columns at the end
    preferred_present = [c for c in spec.column_order if c in out.columns]
    remaining = [c for c in out.columns if c not in preferred_present]
    out = out[preferred_present + remaining]

    # Round configured float columns
    existing_floats = [c for c in spec.float_columns if c in out.columns]
    if existing_floats:
        out[existing_floats] = (
            out[existing_floats]
            .apply(pd.to_numeric, errors="coerce")
            .round(spec.decimals)
        )

    # Stable genomic sort if the required interval columns exist
    chr_col, start_col, end_col = spec.sort_keys
    if all(c in out.columns for c in (chr_col, start_col, end_col)):
        out = _sort_by_chr_interval(out, chr_col, start_col, end_col)

    return out


def _sort_by_chr_interval(
    df: pd.DataFrame,
    chr_col: str,
    start_col: str,
    end_col: str,
    tmp_col: str = "chr_sort",
) -> pd.DataFrame:
    """Stable-sort rows by chromosome, start, and end."""
    out = df.copy()
    out[tmp_col] = out[chr_col].map(chrom_sort_key)
    out = out.sort_values(by=[tmp_col, start_col, end_col], kind="stable").drop(
        columns=[tmp_col]
    )
    return out


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

    if chrom.isdigit():
        expected_cn = 2
    elif chrom.upper() == "X":
        expected_cn = 2 if sex == Gender.FEMALE else 1
    elif chrom.upper() == "Y":
        expected_cn = 1 if sex == Gender.MALE else 0
    else:
        expected_cn = 2

    if cn < expected_cn:
        return "DELETION"
    if cn > expected_cn:
        return "AMPLIFICATION"
    return "NEUTRAL"


def add_sex_aware_cnv_calls_from_total_cn(
    df: pd.DataFrame,
    sex: Gender,
    chr_col: str = CHR,
) -> pd.DataFrame:
    """
    Add CNV gain/loss classifications derived from absolute copy number.

    This function computes CNV calls separately for CNVkit and PureCN
    segments using `classify_cnv_from_total_cn_sex_aware`, which interprets
    copy number relative to the expected baseline for each chromosome and sex.

    Baseline expectations:
        autosomes (1–22): 2 copies
        chromosome X:     2 copies (female) / 1 copy (male)
        chromosome Y:     0 copies (female) / 1 copy (male)

    If the required copy-number column exists, the function adds a new column
    containing the derived CNV call for that caller.

    Parameters
    ----------
    df
        Input dataframe containing segment rows.
    sex
        Sample sex used for sex-chromosome baseline interpretation.


    Returns
    -------
    pd.DataFrame
        Copy of `df` with CNV classification columns added.
    """
    out = df.copy()

    # Column names
    cnvkit_total_cn_col: str = "cnvkit_seg_cn"
    purecn_total_cn_col: str = "purecn_C"
    out_cnvkit_col: str = "cnvkit_cnv_call"
    out_purecn_col: str = "purecn_cnv_call"

    # Create output columns so downstream code can rely on their existence
    out[out_cnvkit_col] = pd.NA
    out[out_purecn_col] = pd.NA

    out[out_cnvkit_col] = out.apply(
        lambda r: classify_cnv_from_total_cn_sex_aware(
            r.get(cnvkit_total_cn_col), r[chr_col], sex
        ),
        axis=1,
    )

    out[out_purecn_col] = out.apply(
        lambda r: classify_cnv_from_total_cn_sex_aware(
            r.get(purecn_total_cn_col), r[chr_col], sex
        ),
        axis=1,
    )

    return out


############################
# GENE REGION LEVEL
############################


def annotate_regions_with_overlapping_segments(
    regions_df: pd.DataFrame,
    segs_df: pd.DataFrame,
    field_map: Mapping[str, str] | None = None,
    seg_start_col: str = "start",
    seg_end_col: str = "end",
) -> pd.DataFrame:
    """
    Generic annotator: for each region, pick the segment with max bp overlap and copy fields.

    field_map maps segs_df column -> regions_df output column name.
    seg_*_col allow segment coordinate columns to be prefixed (e.g. cnvkit_seg_start).
    """
    out = regions_df.copy()

    # Stable column names
    region_chr_col: str = CHR
    region_start_col: str = "region_start"
    region_end_col: str = "region_end"
    seg_chr_col: str = "chr"

    # Keep only fields present in segs_df
    field_map = {src: dst for src, dst in field_map.items() if src in segs_df.columns}

    # Per-chromosome lookup for repeated access inside the region loop
    segs_by_chr: dict[str, pd.DataFrame] = {
        str(ch): df for ch, df in segs_df.groupby(seg_chr_col, sort=False)
    }

    # Ensure destination columns exist
    for src, dst in field_map.items():
        if dst not in out.columns:
            out[dst] = pd.Series(pd.NA, index=out.index, dtype=segs_df[src].dtype)

    # Main loop
    for i, row in out.iterrows():
        ch = str(row.get(region_chr_col, ""))
        segs_chr = segs_by_chr.get(ch)
        if segs_chr is None or segs_chr.empty:
            continue

        picked = _pick_best_overlapping_segment(
            segs_chr,
            row.get(region_start_col, np.nan),
            row.get(region_end_col, np.nan),
            seg_start_col=seg_start_col,
            seg_end_col=seg_end_col,
        )
        if picked is None:
            continue

        for src, dst in field_map.items():
            val = picked.get(src, pd.NA)
            out.at[i, dst] = val

    return out


def _pick_best_overlapping_segment(
    segs_chr: pd.DataFrame,
    region_start: float,
    region_end: float,
    seg_start_col: str = "start",
    seg_end_col: str = "end",
) -> pd.Series | None:
    """
    Return the segment with the largest base-pair overlap with a region.

    The function scans all segments on a chromosome and selects the one whose
    genomic interval overlaps the region (region_start, region_end) by the
    greatest number of base pairs.

    Overlap length is computed as:

        overlap = min(segment_end, region_end) - max(segment_start, region_start)

    Only segments with positive overlap are considered.

    Parameters
    ----------
    segs_chr
        DataFrame containing segment intervals for a single chromosome.
    region_start
        Start coordinate of the region of interest.
    region_end
        End coordinate of the region of interest.
    seg_start_col
        Column name containing segment start coordinates.
    seg_end_col
        Column name containing segment end coordinates.

    Returns
    -------
    pd.Series | None
        The row corresponding to the segment with the largest overlap with the
        region. Returns None if no segments overlap the region.
    """

    if segs_chr is None or segs_chr.empty:
        return None

    overlaps = (segs_chr[seg_end_col] > region_start) & (
        segs_chr[seg_start_col] < region_end
    )
    if not overlaps.any():
        return None

    overlapping_segs = segs_chr.loc[overlaps]
    overlap_start = np.maximum(overlapping_segs[seg_start_col].to_numpy(), region_start)
    overlap_end = np.minimum(overlapping_segs[seg_end_col].to_numpy(), region_end)
    overlap_lengths = overlap_end - overlap_start
    best_idx = int(overlap_lengths.argmax())
    return overlapping_segs.iloc[best_idx]


def aggregate_gene_bins(bins: pd.DataFrame) -> pd.DataFrame:
    """Aggregate bin-level rows to one row per gene-region."""
    return (
        bins.groupby([CHR, GENE], as_index=False)
        .agg(
            region_start=("start", "min"),
            region_end=("end", "max"),
            **{"n.targets": ("log2", "count")},
            mean_log2=("log2", "mean"),
            min_log2=("log2", "min"),
            max_log2=("log2", "max"),
        )
        .copy()
    )


def build_generegion_table(
    cnr_df: pd.DataFrame,
    cns_df: pd.DataFrame,
    sex: Gender,
    pon_df: pd.DataFrame | None = None,
    cancer_genes: set[str] | None = None,
    loh_segments_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """
    Build per-gene table always.
    If pon_df exists -> subdivided into regions based on PON
    Else -> one-row-per-gene 'genelevel'.
    """

    cancer_genes = cancer_genes or set()

    # If no PON exists, just merge CNR BINS per gene level
    if pon_df is None:
        regions_df = aggregate_gene_bins(cnr_df)
    else:
        regions_df = create_generegions(cnr_df=cnr_df, pon_df=pon_df)

    # Annotate with CNVkit segments
    regions_df = annotate_regions_with_overlapping_segments(
        regions_df,
        cns_df,
        field_map={c: c for c in cns_df.columns if c.startswith("cnvkit_")},
        seg_start_col="cnvkit_seg_start",
        seg_end_col="cnvkit_seg_end",
    )

    # Attach PureCN LOHregions (optional)
    if loh_segments_df is not None and not loh_segments_df.empty:
        regions_df = annotate_regions_with_overlapping_segments(
            regions_df,
            loh_segments_df,
            field_map={
                c: c for c in loh_segments_df.columns if c.startswith("purecn_")
            },
            seg_start_col=f"purecn_seg_start",
            seg_end_col=f"purecn_seg_end",
        )

    regions_df["is_cancer_gene"] = regions_df[GENE].isin(cancer_genes)

    # CNV calls
    regions_df = add_sex_aware_cnv_calls_from_total_cn(regions_df, sex=sex)

    # Remove columns
    regions_df = regions_df.drop(
        columns=[
            "cnvkit_seg_start",
            "cnvkit_seg_end",
            "cnvkit_seg_log2",
            "purecn_seg_start",
            "purecn_seg_end",
            "purecn_seg_mean_log2",
            "purecn_seg_mean_log2",
            "purecn_num_snps",
            "purecn_maf_observed",
            "purecn_M_flagged",
            "region_id",
            "pon_region_direction",
            "cnvkit_seg_depth",
        ],
        errors="ignore",
    )

    return finalize_table(regions_df, GENE_TABLE_SPEC)


############################
# SEGMENT LEVEL
############################


def add_overlapping_genes_from_bins(
    seg_df: pd.DataFrame,
    cnr_df: pd.DataFrame,
    *,
    genes_col: str = GENE,
    out_targets_col: str = "n.targets",
    min_targets: int = 2,
    cancer_genes: set[str] | None = None,
) -> pd.DataFrame:
    """
    For each segment row (chr/start/end), find overlapping CNVkit bins in cnr_df and:
      1) add comma-separated gene list in `genes_col`
      2) add `out_targets_col` = number of UNIQUE overlapping bins

    Notes
    -----
    - `cnr_df` is typically exploded by gene.symbol, so `out_targets_col` counts
      unique bins by (chr, start, end) to avoid overcounting.
    - The gene list excludes "backbone".
    - Only genes with at least `min_targets` overlapping bins are included.
    - If `cancer_genes` is provided, only those genes are shown in the gene list.
    """
    out = seg_df.copy()

    # Prepare output columns
    out[genes_col] = ""
    out[out_targets_col] = 0

    # Work per chromosome
    for chrom, seg_g in out.groupby(CHR, sort=False):
        bins_g = cnr_df[cnr_df[CHR] == chrom].copy()

        bins_g = bins_g.sort_values("start", kind="stable").reset_index(drop=True)

        b_start = bins_g["start"].to_numpy(dtype=np.int64, copy=False)
        b_end = bins_g["end"].to_numpy(dtype=np.int64, copy=False)
        b_gene = bins_g[GENE].to_numpy(dtype=object, copy=False)

        seg_idx = seg_g.index.to_numpy()
        s_start = seg_g["start"].to_numpy(dtype=np.int64, copy=False)
        s_end = seg_g["end"].to_numpy(dtype=np.int64, copy=False)

        for i, row_idx in enumerate(seg_idx):
            lo = s_start[i]
            hi = s_end[i]
            if hi <= lo:
                continue

            cut = np.searchsorted(b_start, hi, side="left")
            if cut == 0:
                continue

            mask = b_end[:cut] > lo
            if not np.any(mask):
                continue

            # ------------------------------------------------------
            # (A) n.targets = all UNIQUE overlapping bins in segment
            # ------------------------------------------------------
            starts = b_start[:cut][mask]
            ends = b_end[:cut][mask]
            uniq_bins = np.unique(np.stack([starts, ends], axis=1), axis=0)
            out.at[row_idx, out_targets_col] = int(uniq_bins.shape[0])

            # ------------------------------------------------------
            # (B) gene list = overlapping genes, excluding backbone,
            # ------------------------------------------------------
            genes_overlap = (
                pd.Series(b_gene[:cut][mask], dtype="string").fillna("").astype(str)
            )
            genes_overlap = genes_overlap[genes_overlap.str.strip().ne("")]

            genes_overlap = genes_overlap[
                ~genes_overlap.str.strip().str.lower().eq("backbone")
            ]

            # ------------------------------------------------------
            # (C) If exome, also include only cancer genes
            # ------------------------------------------------------
            if cancer_genes:
                genes_overlap = genes_overlap[genes_overlap.isin(cancer_genes)]

            if genes_overlap.empty:
                continue

            # ------------------------------------------------------
            # (D) Keep only genes above min number of targets
            # ------------------------------------------------------
            gene_bin_counts = genes_overlap.value_counts(dropna=True)
            gene_bin_counts = gene_bin_counts[gene_bin_counts >= min_targets]
            if gene_bin_counts.empty:
                continue

            gene_list = sorted(gene_bin_counts.index.astype(str).tolist())
            out.at[row_idx, genes_col] = ",".join(gene_list)

    return out


def annotate_segments_with_cytoband(
    df_segments: pd.DataFrame,
    cyto: pd.DataFrame,
    *,
    seg_chr_col: str = CHR,
    seg_start_col: str = "start",
    seg_end_col: str = "end",
    cyto_chr_col: str = CHR,
    cyto_start_col: str = "start_int",
    cyto_end_col: str = "end_int",
    cyto_name_col: str = "name",
) -> pd.DataFrame:
    """
    Add 'cytoband' based on overlap between segment (chr,start,end) and cytobands.
    """
    df = df_segments.copy()
    df["cytoband"] = pd.Series(pd.NA, index=df.index, dtype="string")

    # make sure types are sane
    df[seg_chr_col] = df[seg_chr_col].astype("string")
    df[seg_start_col] = pd.to_numeric(df[seg_start_col], errors="coerce").astype(
        "Int64"
    )
    df[seg_end_col] = pd.to_numeric(df[seg_end_col], errors="coerce").astype("Int64")

    cy = cyto.copy()
    cy[cyto_chr_col] = cy[cyto_chr_col].astype("string")
    cy[cyto_start_col] = pd.to_numeric(cy[cyto_start_col], errors="coerce").astype(
        "Int64"
    )
    cy[cyto_end_col] = pd.to_numeric(cy[cyto_end_col], errors="coerce").astype("Int64")

    for chrom, seg_chr in df.groupby(seg_chr_col, sort=False):
        cyto_chr = cy[cy[cyto_chr_col] == chrom]
        if cyto_chr.empty:
            continue

        cyto_chr = cyto_chr.sort_values(cyto_start_col, kind="stable")
        band_starts = cyto_chr[cyto_start_col].to_numpy()
        band_ends = cyto_chr[cyto_end_col].to_numpy()
        band_names = cyto_chr[cyto_name_col].to_numpy()

        # iterate ONLY rows on this chromosome
        for idx, row in seg_chr.iterrows():
            if pd.isna(row[seg_start_col]) or pd.isna(row[seg_end_col]):
                continue

            s_start = int(row[seg_start_col])
            s_end = int(row[seg_end_col])

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


def build_segment_table(
    cnr_df: pd.DataFrame,
    cns_df: pd.DataFrame,
    cytoband_df: pd.DataFrame,
    sex: Gender,
    cancer_genes: set[str],
    is_exome: bool,
    loh_segments_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """
    Build a unified segment table by stacking CNVkit and PureCN segments.

    Output core columns:
      - chr, start, end, caller

    CNVkit-specific (if present in cns_df):
      - cnvkit_seg_cn, cnvkit_seg_cn1, cnvkit_seg_cn2, cnvkit_seg_depth
      - cnvkit_seg_log2, cnvkit_seg_raw_log2

    PureCN-specific (if present in loh_regions_df):
      - purecn_C, purecn_M, purecn_M_flagged, purecn_seg_mean_log2, purecn_num_snps
      - purecn_maf_observed, purecn_type
      - purecn_is_loh (derived: purecn_type contains 'LOH')

    Notes:
      - Duplicated chr/start/end rows are allowed (CNVkit and PureCN can overlap).
      - This function *stacks* callers (long format). It does not attempt interval intersection.
    """
    out_parts: list[pd.DataFrame] = []

    # ---- CNVkit ----
    cnv = cns_df.copy()
    # normalize coords -> shared names
    cnv = cnv.rename(columns={"cnvkit_seg_start": "start"})
    cnv = cnv.rename(columns={"cnvkit_seg_end": "end"})
    cnv = cnv.rename(columns={"cnvkit_seg_log2": "cnvkit_adjusted_log2"})
    cnv = cnv.rename(columns={"cnvkit_seg_raw_log2": "log2"})
    cnv = cnv.rename(columns={"cnvkit_seg_baf": "baf_maf"})

    cnv["caller"] = "CNVkit"
    out_parts.append(cnv)

    # ---- PureCN ----
    if loh_segments_df is not None and not loh_segments_df.empty:
        pc = loh_segments_df.copy()

        pc = pc.rename(columns={"purecn_seg_start": "start"})
        pc = pc.rename(columns={"purecn_seg_end": "end"})
        pc = pc.rename(columns={"purecn_seg_mean_log2": "log2"})
        pc = pc.rename(columns={"purecn_maf_observed": "baf_maf"})

        pc["caller"] = "PureCN"

        out_parts.append(pc)

    segments = pd.concat(out_parts, axis=0, ignore_index=True, sort=False)
    segments = segments.sort_values(
        [CHR, "start", "end", "caller"], kind="stable"
    ).reset_index(drop=True)

    # Add overlapping genes (limit to cancer genes only if provided)
    segments = add_overlapping_genes_from_bins(
        segments,
        cnr_df,
        genes_col=GENE,
        min_targets=2,
        cancer_genes=cancer_genes if is_exome else None,
    )
    segments = add_sex_aware_cnv_calls_from_total_cn(
        segments, sex=sex
    )  # makes cnvkit_cnv_call and purecn_cnv_call

    # Then unify
    segments["cnv_call"] = pd.NA
    mask_cnv = segments["caller"].astype("string").str.upper().eq("CNVKIT")
    mask_pc = segments["caller"].astype("string").str.upper().eq("PURECN")

    segments.loc[mask_cnv, "cnv_call"] = segments.loc[mask_cnv, "cnvkit_cnv_call"]
    segments.loc[mask_pc, "cnv_call"] = segments.loc[mask_pc, "purecn_cnv_call"]

    segments = annotate_segments_with_cytoband(segments, cytoband_df)
    segments["segment_size"] = round((segments["end"] - segments["start"]) / 1000, 2)

    segments = segments.drop(columns=["cnvkit_cnv_call", "purecn_cnv_call"])
    return finalize_table(segments, SEGMENT_TABLE_SPEC)
