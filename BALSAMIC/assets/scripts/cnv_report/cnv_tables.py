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
    GeneRegionConfig,
)
from cnv_report_utils import chrom_sort_key


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
    chr_col: str = "chr",
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


def _pon_abs_z(
    pon_region_log2_difference: float,
    spread: float,
    n_targets: float,
    *,
    min_n: int,
) -> float:
    """
    Compute a PON-relative absolute z-like score for a region.

    The score is defined as:

        abs(pon_region_log2_difference) / (spread / sqrt(n_targets))

    where:
      - pon_region_log2_difference is the region's mean log2 deviation from
        the PON baseline
      - spread is the expected PON variation
      - n_targets is the number of bins contributing to the region

    Returns NaN when the input values are missing, the spread is not positive,
    or the region is too small to score.
    """
    if pd.isna(pon_region_log2_difference) or pd.isna(spread) or spread <= 0:
        return np.nan

    if n_targets < min_n:
        return np.nan

    return abs(pon_region_log2_difference) / (spread / np.sqrt(float(n_targets)))


def _pon_direction(log2difference: float) -> str:
    """Map log2difference sign to 'gain'/'loss'/'neutral' (or '' when missing)."""
    if pd.isna(log2difference):
        return ""
    if log2difference > 0:
        return "gain"
    if log2difference < 0:
        return "loss"
    return "neutral"


def _pon_signal(z: float, *, noise_lt: float, borderline_lt: float) -> str:
    """Map z-score to 'noise'/'borderline'/'strong' (or '' when missing)."""
    if pd.isna(z):
        return ""
    if z < noise_lt:
        return "noise"
    if z < borderline_lt:
        return "borderline"
    return "strong"


def _pon_cnv_call_from_pon_log2_difference(
    *,
    is_strong: bool,
    pon_region_log2_difference: float,
) -> str:
    """Return GAIN/LOSS/neutral/weak depending on deviation from PON baseline."""
    if not is_strong or pd.isna(pon_region_log2_difference):
        return (GeneRegionConfig.pon_neutral_call,)
    if pon_region_log2_difference > GeneRegionConfig.pon_gain_gt:
        return GeneRegionConfig.pon_gain_call
    if pon_region_log2_difference < GeneRegionConfig.pon_loss_lt:
        return GeneRegionConfig.pon_loss_call
    return (GeneRegionConfig.pon_neutral_call,)


def annotate_regions_with_overlapping_segments(
    regions_df: pd.DataFrame,
    segs_df: pd.DataFrame,
    field_map: Mapping[str, str] | None = None,
    *,
    region_chr_col: str = "chr",
    region_start_col: str = "region_start",
    region_end_col: str = "region_end",
    seg_chr_col: str = "chr",
    seg_start_col: str = "start",
    seg_end_col: str = "end",
) -> pd.DataFrame:
    """
    Generic annotator: for each region, pick the segment with max bp overlap and copy fields.

    field_map maps segs_df column -> regions_df output column name.
    seg_*_col allow segment coordinate columns to be prefixed (e.g. cnvkit_seg_start).
    """
    out = regions_df.copy()

    # Keep only fields present in segs_df
    field_map = {src: dst for src, dst in field_map.items() if src in segs_df.columns}

    # Per-chrom map for speed
    segs_by_chr: dict[str, pd.DataFrame] = {
        str(ch): df.reset_index(drop=True)
        for ch, df in segs_df.groupby(seg_chr_col, sort=False)
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
    *,
    seg_start_col: str = "start",
    seg_end_col: str = "end",
) -> pd.Series | None:
    """Return the segment row (Series) with max bp overlap vs region, or None."""
    if segs_chr is None or segs_chr.empty:
        return None
    if pd.isna(region_start) or pd.isna(region_end):
        return None

    mask = (segs_chr[seg_end_col] > region_start) & (
        segs_chr[seg_start_col] < region_end
    )
    if not mask.any():
        return None

    ssub = segs_chr.loc[mask]
    ov0 = np.maximum(ssub[seg_start_col].to_numpy(), region_start)
    ov1 = np.minimum(ssub[seg_end_col].to_numpy(), region_end)
    best = int((ov1 - ov0).argmax())
    return ssub.iloc[best]


def annotate_regions_with_cnvkit_segments(
    regions_df: pd.DataFrame,
    cns_df: pd.DataFrame,
    prefix: str = "cnvkit_",
) -> pd.DataFrame:
    seg_start = f"{prefix}seg_start"
    seg_end = f"{prefix}seg_end"

    field_map = {c: c for c in cns_df.columns if c.startswith(prefix)}

    return annotate_regions_with_overlapping_segments(
        regions_df,
        cns_df,
        field_map=field_map,
        seg_start_col=seg_start,
        seg_end_col=seg_end,
    )


def annotate_regions_with_purecn_lohregions(
    regions_df: pd.DataFrame,
    lohregions_df: pd.DataFrame,
    prefix: str = "purecn_",
) -> pd.DataFrame:
    seg_start = f"{prefix}seg_start"
    seg_end = f"{prefix}seg_end"

    field_map = {c: c for c in lohregions_df.columns if c.startswith(prefix)}

    return annotate_regions_with_overlapping_segments(
        regions_df,
        lohregions_df,
        field_map=field_map,
        seg_start_col=seg_start,
        seg_end_col=seg_end,
    )


def merge_cnr_with_pon(
    cnr_df: pd.DataFrame,
    pon_df: pd.DataFrame | None,
    *,
    key_cols: tuple[str, str, str, str] = ("chr", "start", "end", "gene.symbol"),
    pon_cols: list[str] | None = None,
) -> pd.DataFrame:
    """
    Merge exploded CNR bins with exploded PON bins on (chr,start,end,gene.symbol).

    - Drops backbone rows from both.
    - De-duplicates PON on the merge key to avoid many-to-many inflation.
    - Keeps CNR rows even if PON is missing (left join).

    If pon_cols is given, only those PON columns (plus the key) are merged in.
    """
    cnr = cnr_df.copy()
    pon = pon_df.copy() if pon_df is not None else None

    g_cnr = cnr["gene.symbol"].astype("string").str.strip()
    cnr = cnr.loc[g_cnr.ne("backbone")].copy()

    if pon is None:
        return cnr

    g_pon = pon["gene.symbol"].astype("string").str.strip()
    pon = pon.loc[g_pon.ne("backbone")].copy()

    # Choose which PON columns to bring (prevents collisions like log2/depth/etc.)
    key = list(key_cols)
    if pon_cols is not None:
        keep = [c for c in (key + pon_cols) if c in pon.columns]
        pon = pon[keep].copy()

    # Ensure PON is unique per key (critical for stable merge)
    pon = pon.drop_duplicates(subset=key)

    return cnr.merge(
        pon,
        how="left",
        on=key,
        validate="many_to_one",  # cnr may repeat; pon must be unique per key
    )


def _aggregate_gene_regions(bins: pd.DataFrame) -> pd.DataFrame:
    """Aggregate bin-level rows to one row per gene-region."""
    return (
        bins.groupby(["chr", "gene.symbol", "region_id"], as_index=False)
        .agg(
            region_start=("start", "min"),
            region_end=("end", "max"),
            n_targets=("log2", "count"),
            mean_log2=("log2", "mean"),
            min_log2=("log2", "min"),
            max_log2=("log2", "max"),
            pon_mean_log2=("pon_log2", "mean"),
            pon_mean_spread=("pon_spread", "mean"),
        )
        .copy()
    )


def _prepare_generegion_bins(
    cnr_df: pd.DataFrame,
    pon_df: pd.DataFrame | None,
    *,
    key_cols: tuple[str, str, str, str] = ("chr", "start", "end", "gene.symbol"),
) -> pd.DataFrame:
    """
    Prepare per-bin input table for gene-region construction.

    If PON data is available, merge sample bins with matching PON bins on the
    genomic bin key. If no PON is available, return the sample bins unchanged
    but add empty PON columns so downstream code can use a uniform schema.

    Parameters
    ----------
    cnr_df
        Expanded sample bin table.
    pon_df
        Expanded PON bin table, or None.
    key_cols
        Columns used to match sample bins to PON bins.

    Returns
    -------
    pd.DataFrame
        Sorted bin table with guaranteed columns:
            chr, gene.symbol, start, end, log2, pon_log2, pon_spread
    """

    pon_cols = ["pon_log2", "pon_spread"]
    has_pon = pon_df is not None and not pon_df.empty

    if has_pon:
        bins = merge_cnr_with_pon(
            cnr_df=cnr_df,
            pon_df=pon_df,
            key_cols=key_cols,
            pon_cols=pon_cols,
        )
    else:
        bins = cnr_df.copy()

        # Ensure downstream code can rely on these columns existing
        if "pon_log2" not in bins.columns:
            bins["pon_log2"] = np.nan
        if "pon_spread" not in bins.columns:
            bins["pon_spread"] = np.nan

    bins = bins.sort_values(
        ["chr", "gene.symbol", "start"],
        kind="stable",
    ).reset_index(drop=True)

    return bins


def _build_genelevel_regions_without_pon(bins: pd.DataFrame) -> pd.DataFrame:
    """
    Build one gene-level region per gene when no PON is available.
    """
    gene_bins = bins.copy()
    gene_bins["region_id"] = "genelevel"

    regions_df = _aggregate_gene_regions(gene_bins)

    regions_df["n.targets"] = regions_df["n_targets"]

    # No PON-derived region calling is possible in this mode
    regions_df["pon_region_log2_difference"] = np.nan
    regions_df["pon_region_z"] = np.nan
    regions_df["pon_region_signal"] = ""
    regions_df["pon_region_indication"] = ""

    return regions_df


def _assign_initial_gene_regions(
    bins: pd.DataFrame,
    *,
    config: GeneRegionConfig,
) -> pd.DataFrame:
    """
    Assign provisional gene-region IDs based on runs of bins that deviate
    consistently from the PON baseline.

    For each gene:
      1. Compute per-bin log2difference = log2 - pon_log2
      2. Convert to a per-bin z-like score using pon_spread
      3. Smooth the z-scores across neighboring bins
      4. Find consecutive runs of bins with:
           - sufficiently large absolute z
           - consistent sign (gain-like or loss-like)
      5. Keep only runs that are large enough and strong enough overall

    Bins that do not belong to any qualifying run remain "no_region".
    """
    # Work on a copy so the input dataframe is not modified in-place
    out = bins.copy()

    # Default label for bins that do not end up in any provisional region
    out["region_id"] = "no_region"

    # Global counter so each accepted run gets a unique provisional ID
    next_region_id = 0

    def _run_abs_z(df_run: pd.DataFrame) -> float:
        """
        Compute a run-level absolute z-score.

        This summarizes the full run, rather than judging bins one by one.
        A run is strong when:
          - its mean deviation from PON is large
          - compared with its expected noise
          - taking run length into account
        """
        # Per-bin deviation from the PON baseline
        log2difference = (
            df_run["log2"].to_numpy() - df_run["pon_log2"].fillna(0.0).to_numpy()
        )
        if log2difference.size == 0:
            return np.nan

        # Average deviation across the whole run
        mean_log2difference = float(np.nanmean(log2difference))
        if not np.isfinite(mean_log2difference):
            return np.nan

        # Use the mean PON spread as a rough noise estimate for the run
        spread = df_run["pon_spread"].fillna(0.0).to_numpy()

        # Avoid divide-by-zero for bins with zero or missing spread
        spread_safe = np.where(spread <= 0, 1e-3, spread)

        mean_spread = float(np.nanmean(spread_safe)) if spread_safe.size else np.nan
        if not np.isfinite(mean_spread) or mean_spread <= 0:
            return np.nan

        # Standard error of the mean log2difference across the run
        sigma_log2difference = mean_spread / np.sqrt(float(len(log2difference)))
        if sigma_log2difference <= 0:
            return np.nan

        # Absolute run-level z-score
        return abs(mean_log2difference) / sigma_log2difference

    def _label_gene_runs(gene_bins: pd.DataFrame) -> None:
        """
        Detect and label provisional regions within one gene.

        This function scans bins in genomic order and accumulates runs of bins
        whose smoothed z-scores are both:
          - strong enough in magnitude
          - consistent in direction (positive or negative)
        """
        nonlocal next_region_id

        # Skip tiny genes: too few bins to robustly call regions
        if gene_bins.shape[0] < config.min_gene_targets:
            return

        # Skip genes where there is no usable PON spread at all
        if gene_bins["pon_spread"].isna().all():
            return

        # Make sure bins are processed in genomic order
        gene_bins = gene_bins.sort_values("start", kind="stable")

        # Preserve original dataframe indices so accepted runs can be written
        # back into the full output table
        bin_indices = gene_bins.index.to_numpy()

        # Per-bin log2 difference relative to PON
        log2difference = (
            gene_bins["log2"].to_numpy() - gene_bins["pon_log2"].fillna(0.0).to_numpy()
        )

        # Per-bin spread estimate from PON
        spread = gene_bins["pon_spread"].fillna(0.0).to_numpy()

        # Avoid division by zero when spread is missing/zero
        spread_safe = np.where(spread <= 0, 1e-3, spread)

        # Per-bin z-like deviation score
        z_raw = log2difference / spread_safe

        # Smooth z across neighboring bins to reduce single-bin noise spikes
        z_smooth = (
            pd.Series(z_raw, index=gene_bins.index)
            .rolling(
                window=config.smooth_window,
                center=True,
                min_periods=1,
            )
            .median()
            .to_numpy()
        )

        # Current candidate run of bins
        current_run: list[int] = []

        # Sign of current run:
        #   +1 for gain-like bins
        #   -1 for loss-like bins
        current_sign: int | None = None

        def _finalize_run(run_indices: list[int]) -> None:
            """
            Decide whether the current candidate run is strong enough to keep.

            A run is accepted only if:
              - it contains enough bins
              - its aggregated run-level z-score exceeds threshold
            """
            nonlocal next_region_id

            # Ignore short runs
            if len(run_indices) < config.min_run_bins:
                return

            # Re-evaluate the whole run as one region
            df_run = out.loc[run_indices]
            run_z = _run_abs_z(df_run)

            # Ignore weak or invalid runs
            if not np.isfinite(run_z) or run_z < config.z_run_thresh:
                return

            # Assign a provisional region label to all bins in the run
            region_label = f"generegion_{next_region_id}"
            next_region_id += 1
            out.loc[run_indices, "region_id"] = region_label

        # ------------------------------------------------------------------
        # Main scan across bins:
        # - bins with weak/invalid z break the current run
        # - bins with strong z continue the run if the sign matches
        # - sign changes force the current run to end and a new one to begin
        # ------------------------------------------------------------------
        for bin_index, z_value in zip(bin_indices, z_smooth):
            # Weak or invalid bins cannot belong to a provisional region
            if not np.isfinite(z_value) or abs(z_value) < config.z_bin_thresh:
                if current_run:
                    _finalize_run(current_run)
                    current_run = []
                    current_sign = None
                continue

            # Convert z sign into gain-like (+1) or loss-like (-1)
            z_sign = 1 if z_value > 0 else -1

            # Start a new run, or extend the existing run if sign is consistent
            if current_sign is None or z_sign == current_sign:
                current_run.append(bin_index)
                current_sign = z_sign
            else:
                # Sign flip: finalize old run, then start a new one
                _finalize_run(current_run)
                current_run = [bin_index]
                current_sign = z_sign

        # Final candidate run may still be open when loop ends
        if current_run:
            _finalize_run(current_run)

    # Process each gene independently
    for (_chrom, _gene), gene_bins in out.groupby(["chr", "gene.symbol"], sort=False):
        _label_gene_runs(gene_bins)

    return out


def _merge_adjacent_gene_regions(
    bins: pd.DataFrame,
    *,
    config: GeneRegionConfig,
) -> pd.DataFrame:
    """
    Merge adjacent provisional gene regions within each gene.

    This cleanup step operates after initial region assignment and reduces
    fragmentation caused by short interruptions or very small neighboring runs.

    Merge rules
    -----------
    1. Bridge merge (A-B-C):
       Merge three consecutive runs when:
         - the middle run has at most `config.max_bridge_bins` bins
         - the outer runs have similar mean log2difference
           (difference <= `config.bridge_delta`)

    2. Small-run merge:
       Merge a run with the previous run when:
         - their mean log2difference are similar
           (difference <= `config.merge_delta`)
         - at least one of the two runs has at most `config.small_seg_n` bins

    log2difference is defined as:
        log2 - pon_log2

    Parameters
    ----------
    bins
        Per-bin table with at least:
            chr, gene.symbol, region_id, log2, pon_log2
    config
        Thresholds controlling cleanup and merging.

    Returns
    -------
    pd.DataFrame
        Copy of `bins` with cleaned `region_id` assignments.
    """
    out = bins.copy()

    def _collect_gene_runs(gene_bins: pd.DataFrame) -> tuple[list[dict], np.ndarray]:
        """
        Convert consecutive identical region_id labels into run records.

        Returns
        -------
        runs
            List of dicts with:
              - positions: positional indices within gene_bins
              - indices: original dataframe row indices
              - mean_log2diff: mean(log2 - pon_log2) over the run
        log2difference_all
            Per-bin log2difference array aligned to gene_bins row order.
        """
        log2difference_all = (
            gene_bins["log2"].to_numpy() - gene_bins["pon_log2"].fillna(0.0).to_numpy()
        )
        region_labels = gene_bins["region_id"].tolist()

        runs: list[dict] = []
        current_positions = [0]
        current_label = region_labels[0]

        for pos in range(1, len(region_labels)):
            if region_labels[pos] == current_label:
                current_positions.append(pos)
            else:
                run_indices = gene_bins.index[current_positions].tolist()
                mean_log2difference = (
                    float(np.nanmean(log2difference_all[current_positions]))
                    if current_positions
                    else np.nan
                )
                runs.append(
                    {
                        "positions": current_positions,
                        "indices": run_indices,
                        "mean_log2diff": mean_log2difference,
                    }
                )
                current_label = region_labels[pos]
                current_positions = [pos]

        if current_positions:
            run_indices = gene_bins.index[current_positions].tolist()
            mean_log2difference = (
                float(np.nanmean(log2difference_all[current_positions]))
                if current_positions
                else np.nan
            )
            runs.append(
                {
                    "positions": current_positions,
                    "indices": run_indices,
                    "mean_log2diff": mean_log2difference,
                }
            )

        return runs, log2difference_all

    def _merge_gene_runs(gene_bins: pd.DataFrame) -> None:
        """
        Apply bridge-merge and small-run merge rules to one gene.
        """
        if gene_bins.empty or gene_bins.shape[0] <= 1:
            return

        gene_bins = gene_bins.sort_values("start", kind="stable")
        runs, log2difference_all = _collect_gene_runs(gene_bins)

        if len(runs) <= 1:
            return

        merged_runs: list[dict] = []
        i = 0

        while i < len(runs):
            # --------------------------------------------------------------
            # Bridge merge: A-B-C -> merge if B is short and A/C are similar
            # --------------------------------------------------------------
            if i <= len(runs) - 3:
                run_a = runs[i]
                run_b = runs[i + 1]
                run_c = runs[i + 2]

                if (
                    len(run_b["indices"]) <= config.max_bridge_bins
                    and np.isfinite(run_a["mean_log2diff"])
                    and np.isfinite(run_c["mean_log2diff"])
                    and abs(run_a["mean_log2diff"] - run_c["mean_log2diff"])
                    <= config.bridge_delta
                ):
                    merged_positions = (
                        run_a["positions"] + run_b["positions"] + run_c["positions"]
                    )
                    merged_indices = (
                        run_a["indices"] + run_b["indices"] + run_c["indices"]
                    )
                    merged_mean_log2difference = (
                        float(np.nanmean(log2difference_all[merged_positions]))
                        if merged_positions
                        else np.nan
                    )

                    merged_runs.append(
                        {
                            "positions": merged_positions,
                            "indices": merged_indices,
                            "mean_log2diff": merged_mean_log2difference,
                        }
                    )
                    i += 3
                    continue

            # --------------------------------------------------------------
            # Small-run merge: merge adjacent similar runs if one is small
            # --------------------------------------------------------------
            current_run = runs[i]

            if merged_runs:
                previous_run = merged_runs[-1]

                should_merge = (
                    np.isfinite(current_run["mean_log2diff"])
                    and np.isfinite(previous_run["mean_log2diff"])
                    and abs(
                        current_run["mean_log2diff"] - previous_run["mean_log2diff"]
                    )
                    <= config.merge_delta
                    and (
                        len(current_run["indices"]) <= config.small_segment_max_bins
                        or len(previous_run["indices"]) <= config.small_segment_max_bins
                    )
                )

                if should_merge:
                    merged_positions = (
                        previous_run["positions"] + current_run["positions"]
                    )
                    merged_indices = previous_run["indices"] + current_run["indices"]
                    merged_mean_log2difference = (
                        float(np.nanmean(log2difference_all[merged_positions]))
                        if merged_positions
                        else np.nan
                    )

                    merged_runs[-1] = {
                        "positions": merged_positions,
                        "indices": merged_indices,
                        "mean_log2diff": merged_mean_log2difference,
                    }
                else:
                    merged_runs.append(current_run)
            else:
                merged_runs.append(current_run)

            i += 1

        # Re-label cleaned runs within this gene
        for run_idx, run in enumerate(merged_runs):
            out.loc[run["indices"], "region_id"] = f"generegion_clean_{run_idx}"

    for (_chrom, _gene), gene_bins in out.groupby(["chr", "gene.symbol"], sort=False):
        _merge_gene_runs(gene_bins)

    return out


def _collapse_bins_to_gene_regions(bins: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse region-labeled bins to one row per gene-region.

    Each output row represents a contiguous gene region and summarizes:
      - genomic span
      - number of contributing bins/targets
      - log2 distribution across bins
      - mean PON baseline and spread
    """
    regions_df = _aggregate_gene_regions(bins)

    # Keep compatibility with older downstream code using n.targets
    regions_df["n.targets"] = regions_df["n_targets"]

    return regions_df


def _score_pon_regions(
    regions_df: pd.DataFrame,
    *,
    config: GeneRegionConfig,
) -> pd.DataFrame:
    """
    Compute PON-based log2 difference, z-score, signal class, and CNV indication
    for each gene-region.

    Scoring is based on deviation from the PON baseline:

        pon_region_log2_difference = mean_log2 - pon_mean_log2

    The absolute PON log2 difference is converted into a z-like score using the mean
    PON spread and the number of bins in the region. Regions with too few
    targets are not eligible for PON-based calls.

    Parameters
    ----------
    regions_df
        Gene-region summary table produced by `_collapse_bins_to_gene_regions`.
    config
        Thresholds controlling region scoring and indication calls.

    Returns
    -------
    pd.DataFrame
        Copy of `regions_df` with added columns:
            pon_region_log2_difference
            pon_region_z
            pon_region_direction
            pon_region_signal
            pon_region_indication
    """
    out = regions_df.copy()

    # Mean deviation from the PON baseline
    out["pon_region_log2_difference"] = out["mean_log2"] - out["pon_mean_log2"]

    # Region-level z-score using target count as effective sample size
    out["pon_region_z"] = out.apply(
        lambda row: _pon_abs_z(
            row["pon_region_log2_difference"],
            row["pon_mean_spread"],
            row.get("n_targets", row.get("n.targets", np.nan)),
            min_n=config.min_region_targets_for_call,
        ),
        axis=1,
    )

    # Direction of deviation relative to PON baseline
    out["pon_region_direction"] = out["pon_region_log2_difference"].apply(
        _pon_direction
    )

    # Qualitative signal class from z-score
    out["pon_region_signal"] = out["pon_region_z"].apply(
        lambda z: _pon_signal(
            z,
            noise_lt=config.pon_signal_noise_lt,
            borderline_lt=config.pon_signal_borderline_lt,
        )
    )

    # Hard gate: small regions are not eligible for PON-based interpretation
    too_small = (
        out["n_targets"].fillna(0).astype(int) < config.min_region_targets_for_call
    )

    out.loc[too_small, "pon_region_signal"] = ""
    out.loc[too_small, "pon_region_indication"] = ""

    # Only eligible regions can receive a gain/loss indication
    eligible = ~too_small
    out.loc[eligible, "pon_region_indication"] = out.loc[eligible].apply(
        lambda row: _pon_cnv_call_from_pon_log2_difference(
            is_strong=str(row.get("pon_region_signal", "")).strip().lower() == "strong",
            pon_region_log2_difference=row.get("pon_region_log2_difference", np.nan),
        ),
        axis=1,
    )

    return out


def create_generegions(
    cnr_df: pd.DataFrame,
    pon_df: pd.DataFrame | None = None,
    *,
    config: GeneRegionConfig = GeneRegionConfig(),
) -> pd.DataFrame:
    """
    Build gene-region summary table from CNR bins.

    With PON:
      - subdivide genes into PON-supported regions
      - score each region relative to the PON baseline

    Without PON:
      - return one gene-level row per gene
      - leave PON-derived outputs empty
    """
    bins = _prepare_generegion_bins(cnr_df=cnr_df, pon_df=pon_df)

    if pon_df is None or pon_df.empty:
        return _build_genelevel_regions_without_pon(bins)

    bins = _assign_initial_gene_regions(bins, config=config)
    bins = _merge_adjacent_gene_regions(bins, config=config)

    regions_df = _collapse_bins_to_gene_regions(bins)
    regions_df = _score_pon_regions(regions_df, config=config)
    return regions_df


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
    regions_df = create_generegions(cnr_df=cnr_df, pon_df=pon_df)

    # Annotate with CNVkit segments
    regions_df = annotate_regions_with_cnvkit_segments(regions_df, cns_df)

    # Attach PureCN LOHregions (optional)
    if loh_segments_df is not None and not loh_segments_df.empty:
        regions_df = annotate_regions_with_purecn_lohregions(
            regions_df, loh_segments_df
        )

    regions_df["is_cancer_gene"] = regions_df["gene.symbol"].isin(cancer_genes)

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
            "n_targets",
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
    genes_col: str = "gene.symbol",
    out_targets_col: str = "n.targets",
    min_targets: int = 2,
    cancer_genes: set[str] | None = None,
    drop_genes: set[str] | None = None,
    # which bins count as "targets" for n.targets
    target_drop_genes: set[str] | None = None,
) -> pd.DataFrame:
    """
    For each segment row (chr/start/end), find overlapping CNVkit bins in cnr_df and:
      1) add comma-separated gene list in `genes_col` (genes w/ >= min_targets bins)
      2) add `out_targets_col` = number of UNIQUE overlapping *target* bins

    Notes:
      - cnr_df is typically exploded by gene.symbol, so `out_targets_col` counts unique bins
        by (chr,start,end) to avoid overcounting.
      - `drop_genes` affects the gene list.
      - `target_drop_genes` affects counting targets (defaults to drop_genes if not provided).
    """
    out = seg_df.copy()

    # Normalize segment types
    out["chr"] = out["chr"].astype("string")
    out["start"] = pd.to_numeric(out["start"], errors="coerce").astype("Int64")
    out["end"] = pd.to_numeric(out["end"], errors="coerce").astype("Int64")

    bins = cnr_df.copy()
    bins["chr"] = bins["chr"].astype("string")
    bins["start"] = pd.to_numeric(bins["start"], errors="coerce").astype("Int64")
    bins["end"] = pd.to_numeric(bins["end"], errors="coerce").astype("Int64")
    bins["gene.symbol"] = bins["gene.symbol"].astype("string")

    # Drop unusable rows early
    out = out.dropna(subset=["chr", "start", "end"]).copy()
    bins = bins.dropna(subset=["chr", "start", "end", "gene.symbol"]).copy()

    # normalize drop sets (lowercase)
    drop_set: set[str] = set()
    if drop_genes:
        drop_set = {str(g).strip().lower() for g in drop_genes if str(g).strip()}

    target_drop_set: set[str] = set(drop_set)
    if target_drop_genes is not None:
        target_drop_set = {
            str(g).strip().lower() for g in target_drop_genes if str(g).strip()
        }

    # Prepare output columns
    out[genes_col] = ""
    out[out_targets_col] = 0

    # Work per chromosome
    for chrom, seg_g in out.groupby("chr", sort=False):
        bins_g = bins[bins["chr"] == chrom].copy()
        if bins_g.empty:
            continue

        bins_g = bins_g.sort_values("start", kind="stable").reset_index(drop=True)

        b_start = bins_g["start"].to_numpy(dtype=np.int64, copy=False)
        b_end = bins_g["end"].to_numpy(dtype=np.int64, copy=False)
        b_gene = bins_g["gene.symbol"].to_numpy(dtype=object, copy=False)
        # for unique-bin counting (avoid exploded overcount)
        b_key_start = b_start
        b_key_end = b_end

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

            cand_ends = b_end[:cut]
            mask = cand_ends > lo
            if not np.any(mask):
                continue

            # ----------------------------
            # (A) n.targets = unique bins overlapping, excluding backbone etc
            # ----------------------------
            genes_overlap = (
                pd.Series(b_gene[:cut][mask], dtype="string").fillna("").astype(str)
            )
            if target_drop_set:
                keep_targets = (
                    ~genes_overlap.str.strip().str.lower().isin(target_drop_set)
                )
            else:
                keep_targets = genes_overlap.ne("")

            if keep_targets.any():
                # (A) n.targets = number of UNIQUE bins overlapping this segment (no exclusions)
                starts = b_key_start[:cut][mask]
                ends = b_key_end[:cut][mask]

                # unique (start,end) pairs => unique bins
                uniq = np.unique(np.stack([starts, ends], axis=1), axis=0)
                out.at[row_idx, out_targets_col] = int(uniq.shape[0])

            # ----------------------------
            # (B) gene list (>= min_targets bins per gene) with optional cancer restriction
            # ----------------------------
            genes_for_list = genes_overlap

            if drop_set:
                genes_for_list = genes_for_list[
                    ~genes_for_list.str.strip().str.lower().isin(drop_set)
                ]

            if cancer_genes:
                cancer_set = {str(g).strip() for g in cancer_genes if str(g).strip()}
                if cancer_set:
                    genes_for_list = genes_for_list[genes_for_list.isin(cancer_set)]

            if genes_for_list.empty:
                continue

            vc = genes_for_list.value_counts(dropna=True)
            vc = vc[vc >= min_targets]
            if vc.empty:
                continue

            gene_list = sorted(vc.index.astype(str).tolist())
            out.at[row_idx, genes_col] = ",".join(gene_list)

    return out


def annotate_segments_with_cytoband(
    df_segments: pd.DataFrame,
    cyto: pd.DataFrame,
    *,
    seg_chr_col: str = "chr",
    seg_start_col: str = "start",
    seg_end_col: str = "end",
    cyto_chr_col: str = "chr",
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
        ["chr", "start", "end", "caller"], kind="stable"
    ).reset_index(drop=True)

    # Add overlapping genes (limit to cancer genes only if provided)
    segments = add_overlapping_genes_from_bins(
        segments,
        cnr_df,
        genes_col="gene.symbol",
        min_targets=2,
        cancer_genes=cancer_genes if is_exome else None,
        drop_genes={"backbone"},
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
