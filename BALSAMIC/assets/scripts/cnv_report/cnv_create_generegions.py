from __future__ import annotations

# Third-party
import numpy as np
import pandas as pd
from dataclasses import dataclass

# Local
from cnv_constants import GeneRegionConfig, GENE, CHR


###############################
# ANNOTATE CNR BINS WITH PON INFO
###############################


def annotate_cnr_bins_with_pon(
    cnr_df: pd.DataFrame,
    pon_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Attach selected PON columns to exploded CNR bins using exact bin-key matching.

    Both inputs are expected to already be expanded on `gene.symbol`. Rows are
    matched on `(chr, start, end, gene.symbol)`, and matching PON columns are
    appended to the CNR rows.

    Backbone rows are removed from both inputs before matching. If `pon_df` is
    None, the filtered CNR table is returned unchanged.
    """
    cnr = cnr_df.copy()
    pon = pon_df.copy()

    g_cnr = cnr[GENE]
    cnr = cnr.loc[g_cnr.ne("backbone")].copy()

    g_pon = pon[GENE].astype("string").str.strip()
    pon = pon.loc[g_pon.ne("backbone")].copy()

    merge_columns = [CHR, "start", "end", GENE]
    pon_cols = ["pon_log2", "pon_spread"]

    keep = [c for c in (merge_columns + pon_cols) if c in pon.columns]
    pon = pon[keep].copy()

    return cnr.merge(
        pon,
        how="left",
        on=merge_columns,
        validate="many_to_one",
    )


###############################
# ASSIGN INITIAL GENE REGIONS
###############################


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


def _assign_initial_gene_regions(
    bins: pd.DataFrame,
) -> pd.DataFrame:
    """
    Assign provisional gene-region IDs based on contiguous runs of bins
    that deviate consistently from the PON baseline.

    Strategy
    --------
    For each gene:
      1. Compute a per-bin PON-relative z-like score
      2. Smooth that score across neighboring bins
      3. Scan bins in genomic order and collect consecutive runs where:
           - absolute z is strong enough
           - the sign is consistent across the run
      4. Re-score each candidate run as a whole
      5. Keep only runs that are long enough and strong enough overall

    Bins that do not belong to any accepted run remain labeled "no_region".
    """
    out = bins.copy()
    out["region_id"] = "no_region"

    next_region_id = 0

    # Process one gene at a time so regions are only formed within genes
    for (_chrom, _gene), gene_bins in out.groupby([CHR, GENE], sort=False):
        accepted_runs = _find_gene_runs(gene_bins)

        for run_indices in accepted_runs:
            out.loc[run_indices, "region_id"] = f"generegion_{next_region_id}"
            next_region_id += 1

    return out


def _find_gene_runs(gene_bins: pd.DataFrame) -> list[list[int]]:
    """
    Find accepted candidate runs within a single gene.

    Returns
    -------
    list[list[int]]
        A list of accepted runs, where each run is a list of original
        dataframe indices belonging to that run.
    """
    # Very small genes are skipped entirely because they do not provide
    # enough bins to robustly define a region.
    if len(gene_bins) < GeneRegionConfig.min_gene_targets:
        return []

    # Work in genomic order so neighboring bins are evaluated correctly.
    gene_bins = gene_bins.sort_values("start", kind="stable")

    # -----------------------------------------
    # The raw z values are median-smoothed across neighboring bins
    # -----------------------------------------
    z_smooth = _compute_smoothed_bin_z(gene_bins)

    bin_indices = gene_bins.index.to_numpy()

    accepted_runs: list[list[int]] = []
    current_run: list[int] = []
    current_sign: int | None = None

    # Walk across bins from left to right.
    for bin_index, z_value in zip(bin_indices, z_smooth):
        # Weak bins cannot belong to a candidate run.
        # They terminate any run currently in progress.
        if not abs(z_value) < GeneRegionConfig.z_bin_thresh:
            if _run_passes_thresholds(gene_bins, current_run):
                accepted_runs.append(current_run.copy())
            current_run = []
            current_sign = None
            continue

        # Convert the z-score sign into gain-like (+1) or loss-like (-1).
        z_sign = 1 if z_value > 0 else -1

        # Extend the current run if the direction is consistent.
        if current_sign is None or z_sign == current_sign:
            current_run.append(bin_index)
            current_sign = z_sign
            continue

        # -----------------------------------------
        # Sign has flipped from + to - (gain to loss)
        # Start new run
        # -----------------------------------------
        if _run_passes_thresholds(gene_bins, current_run):
            accepted_runs.append(current_run.copy())

        current_run = [bin_index]
        current_sign = z_sign

    # Final run may still be open after the loop ends.
    if _run_passes_thresholds(gene_bins, current_run):
        accepted_runs.append(current_run.copy())

    return accepted_runs


def _compute_smoothed_bin_z(gene_bins: pd.DataFrame) -> np.ndarray:
    """
    Compute smoothed per-bin PON-relative z-like scores for one gene.

    Per-bin score:
        z = (log2 - pon_log2) / pon_spread

    The raw z values are then median-smoothed across neighboring bins
    to reduce single-bin noise spikes.
    """
    # Per-bin deviation from the PON baseline
    log2difference = (
        gene_bins["log2"].to_numpy() - gene_bins["pon_log2"].fillna(0.0).to_numpy()
    )

    # Expected variability from the PON.
    spread = gene_bins["pon_spread"].to_numpy()

    z_raw = log2difference / spread

    # Smooth across neighboring bins so isolated spikes are less likely
    # to seed false candidate regions.
    z_smooth = (
        pd.Series(z_raw, index=gene_bins.index)
        .rolling(
            window=GeneRegionConfig.smooth_window,
            center=True,
            min_periods=1,
        )
        .median()
        .to_numpy()
    )

    return z_smooth


def _run_passes_thresholds(
    gene_bins: pd.DataFrame,
    run_indices: list[int],
) -> bool:
    """
    Decide whether a candidate run is strong enough to keep.

    A run is accepted only if:
      - it contains enough bins
      - its run-level absolute z-score exceeds the configured threshold
    """
    if len(run_indices) < GeneRegionConfig.min_run_bins:
        return False

    df_run = gene_bins.loc[run_indices]
    run_z = _run_abs_z(df_run)

    return run_z >= GeneRegionConfig.z_run_thresh


def _run_abs_z(df_run: pd.DataFrame) -> float:
    """
    Compute a run-level absolute z-score.

    This evaluates the run as a whole rather than judging each bin separately.

    Formula
    -------
        abs(mean(log2 - pon_log2)) / (mean(pon_spread) / sqrt(n_bins))

    Interpretation
    --------------
    - numerator: average deviation from the PON baseline across the run
    - denominator: expected uncertainty of that average
    - longer runs become more convincing because the standard error shrinks
    """
    # Mean signal across the run
    log2difference = df_run["log2"].to_numpy() - df_run["pon_log2"].to_numpy()
    mean_log2difference = float(np.nanmean(log2difference))

    # Mean expected variability across the run
    spread = df_run["pon_spread"].to_numpy()
    mean_spread = float(np.nanmean(spread))

    run_std_error = mean_spread / np.sqrt(len(df_run))

    return abs(mean_log2difference) / run_std_error


###############################
# MERGE INITIAL GENE REGIONS
###############################


@dataclass
class GeneRun:
    """Contiguous run of bins sharing the same provisional region label."""

    positions: list[int]  # positional indices within gene_bins row order
    indices: list[int]  # original dataframe indices in `out`
    mean_log2diff: float  # mean(log2 - pon_log2) across the run


def _merge_adjacent_gene_regions(
    bins: pd.DataFrame,
) -> pd.DataFrame:
    """
    Merge adjacent provisional gene regions within each gene.

    This is a cleanup step after initial region assignment. It reduces
    fragmentation caused by short interruptions or very small neighboring runs.

    Merge rules
    -----------
    1. Bridge merge (A-B-C):
       Merge three consecutive runs when:
         - the middle run is short
         - the outer runs have similar mean log2difference

    2. Small-run merge:
       Merge two adjacent runs when:
         - their mean log2difference are similar
         - at least one of the runs is small
    """
    out = bins.copy()

    for (_chrom, _gene), gene_bins in out.groupby([CHR, GENE], sort=False):
        gene_bins = gene_bins.sort_values("start", kind="stable")

        # Merge gene runs per gene if appropriate
        merged_runs = _merge_gene_runs_for_one_gene(gene_bins)

        # Re-label merged runs within this gene
        for run_idx, run in enumerate(merged_runs):
            out.loc[run.indices, "region_id"] = f"generegion_clean_{run_idx}"

    return out


def _merge_gene_runs_for_one_gene(gene_bins: pd.DataFrame) -> list[GeneRun]:
    """
    Collect contiguous runs for one gene and merge them according to the
    bridge-merge rule.
    """

    # Runs: list[GeneRun]
    runs, log2difference_all = _collect_gene_runs(gene_bins)

    if len(runs) <= 1:
        # Only one run for the gene available, just return it
        return runs

    merged_runs: list[GeneRun] = []
    i = 0

    while i < len(runs):
        # --------------------------------------------------------------
        # Bridge merge: A-B-C -> merge if B is short and A/C are similar
        # --------------------------------------------------------------
        # Only do this if there are at least 3 runs left
        if i + 2 < len(runs):
            run_a = runs[i]
            run_b = runs[i + 1]
            run_c = runs[i + 2]

            if _can_bridge_merge(run_a, run_b, run_c):
                merged_runs.append(
                    _combine_runs([run_a, run_b, run_c], log2difference_all)
                )
                i += 3
                continue

        current_run = runs[i]

        merged_runs.append(current_run)
        i += 1

    return merged_runs


def _collect_gene_runs(gene_bins: pd.DataFrame) -> tuple[list[GeneRun], np.ndarray]:
    """
    Convert consecutive identical region_id labels into contiguous run records.
    """
    log2difference_all = gene_bins["log2"].to_numpy() - gene_bins["pon_log2"].to_numpy()
    # Example:
    # ["regionA", "regionA", "regionA", "regionA", "no_region", "regionB", "regionB", "regionB", "regionB"]
    bin_region_labels: list = gene_bins["region_id"].tolist()

    runs: list[GeneRun] = []
    n_rows = len(bin_region_labels)

    run_start = 0

    for pos in range(1, n_rows):
        previous_label = bin_region_labels[pos - 1]
        current_label = bin_region_labels[pos]

        # No new generegion yet: keep extending the current run
        if current_label == previous_label:
            continue

        # New generegion encountered: previous run spans run_start .. pos-1
        run_positions = list(range(run_start, pos))
        runs.append(_build_run(gene_bins, run_positions, log2difference_all))

        # New generegion run starts here
        run_start = pos

    # Add the final run
    run_positions = list(range(run_start, n_rows))
    runs.append(_build_run(gene_bins, run_positions, log2difference_all))

    return runs, log2difference_all


def _build_run(
    gene_bins: pd.DataFrame,
    positions: list[int],
    log2difference_all: np.ndarray,
) -> GeneRun:
    """
    Build one GeneRun from a list of positional indices within gene_bins.
    """
    return GeneRun(
        positions=positions.copy(),
        indices=gene_bins.index[positions].tolist(),
        mean_log2diff=float(np.nanmean(log2difference_all[positions])),
    )


def _combine_runs(
    runs: list[GeneRun],
    log2difference_all: np.ndarray,
) -> GeneRun:
    """
    Merge several runs into one combined run and recompute mean log2difference.
    """
    merged_positions: list[int] = []
    merged_indices: list[int] = []

    for run in runs:
        merged_positions.extend(run.positions)
        merged_indices.extend(run.indices)

    return GeneRun(
        positions=merged_positions,
        indices=merged_indices,
        mean_log2diff=float(np.nanmean(log2difference_all[merged_positions])),
    )


def _can_bridge_merge(run_a: GeneRun, run_b: GeneRun, run_c: GeneRun) -> bool:
    """
    Return True if the runs A-B-C should be merged as a bridge pattern.

    Bridge merge occurs when:
      - the middle run (B) is very short
      - the outer runs (A and C) have similar mean log2difference
    """

    # Middle run must be small
    middle_run_is_small = len(run_b.indices) <= GeneRegionConfig.max_bridge_bins

    # Outer runs must have similar signal strength
    outer_runs_are_similar = (
        abs(run_a.mean_log2diff - run_c.mean_log2diff) <= GeneRegionConfig.bridge_delta
    )

    return middle_run_is_small and outer_runs_are_similar


###################################################
# COLLAPSE BINS INTO ASSIGNED GENE REGIONS
###################################################


def _collapse_bins_into_assigned_gene_regions(bins: pd.DataFrame) -> pd.DataFrame:
    """Aggregate bin-level rows to one row per gene-region."""
    return (
        bins.groupby([CHR, GENE, "region_id"], as_index=False)
        .agg(
            region_start=("start", "min"),
            region_end=("end", "max"),
            **{"n.targets": ("log2", "count")},
            mean_log2=("log2", "mean"),
            min_log2=("log2", "min"),
            max_log2=("log2", "max"),
            pon_mean_log2=("pon_log2", "mean"),
            pon_mean_spread=("pon_spread", "mean"),
        )
        .copy()
    )


###################################################
# SCORE GENE REGIONS ( GAIN / LOSS / NEUTRAL )
###################################################


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
    is_strong: bool,
    pon_region_log2_difference: float,
) -> str:
    """Return GAIN/LOSS/neutral/weak depending on deviation from PON baseline."""
    if not is_strong or pd.isna(pon_region_log2_difference):
        return GeneRegionConfig.pon_neutral_call
    if pon_region_log2_difference > GeneRegionConfig.pon_gain_gt:
        return GeneRegionConfig.pon_gain_call
    if pon_region_log2_difference < GeneRegionConfig.pon_loss_lt:
        return GeneRegionConfig.pon_loss_call
    return GeneRegionConfig.pon_neutral_call


def _score_pon_regions(
    regions_df: pd.DataFrame,
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
            row.get("n.targets"),
            min_n=GeneRegionConfig.min_region_targets_for_call,
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
            noise_lt=GeneRegionConfig.pon_signal_noise_lt,
            borderline_lt=GeneRegionConfig.pon_signal_borderline_lt,
        )
    )

    # Hard gate: small regions are not eligible for PON-based interpretation
    too_small = (
        out["n.targets"].fillna(0).astype(int)
        < GeneRegionConfig.min_region_targets_for_call
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


###################################################
# MAIN FUNCTION
###################################################


def create_generegions(
    cnr_df: pd.DataFrame,
    pon_df: pd.DataFrame,
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
    # Annotate CNR bins with information from overlapping PON bins
    bins = annotate_cnr_bins_with_pon(
        cnr_df=cnr_df,
        pon_df=pon_df,
    )

    # Find bins that differ in consistent way from PON log2 and spread:
    # Assign to a gene region
    bins = _assign_initial_gene_regions(bins)

    # Gene-regions may be subdivided by small noisy bins:
    # Merge suitable gene-regions separated only by few bins
    bins = _merge_adjacent_gene_regions(bins)

    # Collapse bins into their assigned gene regions
    regions_df = _collapse_bins_into_assigned_gene_regions(bins)

    # Score the PON regions
    # Strong signal?
    # Gain / Loss / Neutral?
    regions_df = _score_pon_regions(regions_df)

    return regions_df
