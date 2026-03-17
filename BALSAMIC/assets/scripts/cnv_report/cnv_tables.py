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


def reorder_and_sort_table(df: pd.DataFrame, spec: TableSpec) -> pd.DataFrame:
    """
    Finalize a table based on a TableSpec: rename, reorder, round floats, stable-sort by interval.
    """
    out = df.copy()

    # 1) Reorder columns (keep extras at end)
    out = _reorder_columns(out, list(spec.column_order))

    # 2) Round floats
    existing_floats = [c for c in spec.float_columns if c in out.columns]
    if existing_floats:
        out[existing_floats] = (
            out[existing_floats]
            .apply(pd.to_numeric, errors="coerce")
            .round(spec.decimals)
        )

    # 3) Sort
    chr_col, start_col, end_col = spec.sort_keys
    if all(c in out.columns for c in (chr_col, start_col, end_col)):
        out = _stable_sort_by_chr_interval(out, chr_col, start_col, end_col)

    return out


def _stable_sort_by_chr_interval(
    df: pd.DataFrame,
    chr_col: str,
    start_col: str,
    end_col: str,
    *,
    tmp_col: str = "chr_sort",
) -> pd.DataFrame:
    """Stable-sort by (chr, start, end) using `chrom_sort_key`."""
    df = df.copy()
    df[tmp_col] = df[chr_col].map(chrom_sort_key)
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


def _pon_signal(z: float, *, noise_lt: float, borderline_lt: float) -> str:
    """Map z-score to 'noise'/'borderline'/'strong' (or '' when missing)."""
    if pd.isna(z):
        return ""
    if z < noise_lt:
        return "noise"
    if z < borderline_lt:
        return "borderline"
    return "strong"


def _pon_cnv_call_from_effect(
    *,
    is_strong: bool,
    effect_log2: float,
    gain_gt: float,
    loss_lt: float,
    weak_value: str,
    neutral_value: str,
) -> str:
    """Return GAIN/LOSS/neutral/weak depending on deviation from PON baseline."""
    if not is_strong or pd.isna(effect_log2):
        return weak_value
    if effect_log2 > gain_gt:
        return "GAIN"
    if effect_log2 < loss_lt:
        return "LOSS"
    return neutral_value


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
        expected = 2
    elif chrom.upper() == "X":
        expected = 2 if sex == Gender.FEMALE else 1
    elif chrom.upper() == "Y":
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


def add_cnv_calls_wide(
    df: pd.DataFrame,
    *,
    sex: Gender,
    chr_col: str = "chr",
    cnvkit_total_cn_col: str = "cnvkit_seg_cn",
    purecn_total_cn_col: str = "purecn_C",
    out_cnvkit_col: str = "cnvkit_cnv_call",
    out_purecn_col: str = "purecn_cnv_call",
) -> pd.DataFrame:
    out = df.copy()

    # Always create output columns so downstream code can rely on them
    out[out_cnvkit_col] = pd.NA
    out[out_purecn_col] = pd.NA

    if cnvkit_total_cn_col in out.columns:
        out[out_cnvkit_col] = out.apply(
            lambda r: classify_cnv_from_total_cn_sex_aware(
                r[cnvkit_total_cn_col], r[chr_col], sex
            ),
            axis=1,
        )

    if purecn_total_cn_col in out.columns:
        out[out_purecn_col] = out.apply(
            lambda r: classify_cnv_from_total_cn_sex_aware(
                r[purecn_total_cn_col], r[chr_col], sex
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

    In the no-PON case, genes are not subdivided into multiple regions.
    Instead, all bins belonging to the same gene are collapsed into a single
    row spanning the full gene extent.

    The returned table matches the usual gene-region schema as closely as
    possible, but PON-derived columns are left empty or neutral.

    Parameters
    ----------
    bins
        Per-bin table containing at least:
            chr, gene.symbol, start, end, log2, pon_log2, pon_spread

    Returns
    -------
    pd.DataFrame
        One row per gene with columns including:
            chr, gene.symbol, region_id, region_start, region_end,
            n_targets, mean_log2, min_log2, max_log2,
            pon_mean_log2, pon_mean_spread,
            n.targets,
            pon_region_effect, pon_region_z,
            pon_region_signal, pon_region_indication
    """
    gene_bins = bins.copy()
    gene_bins["region_id"] = "genelevel"

    agg_dict: dict[str, list[str]] = {
        "start": ["min"],
        "end": ["max"],
        "log2": ["count", "mean", "min", "max"],
        "pon_log2": ["mean"],
        "pon_spread": ["mean"],
    }

    regions_df = (
        gene_bins.groupby(["chr", "gene.symbol", "region_id"], as_index=False)
        .agg(agg_dict)
        .copy()
    )

    regions_df = _flatten_agg_columns(regions_df).rename(
        columns={
            "start_min": "region_start",
            "end_max": "region_end",
            "log2_count": "n_targets",
            "log2_mean": "mean_log2",
            "log2_min": "min_log2",
            "log2_max": "max_log2",
            "pon_log2_mean": "pon_mean_log2",
            "pon_spread_mean": "pon_mean_spread",
        }
    )

    # Keep both naming conventions for compatibility with downstream code
    regions_df["n.targets"] = regions_df["n_targets"]

    # No PON-derived region calling is possible in this mode
    regions_df["pon_region_effect"] = np.nan
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
    Assign provisional gene-region IDs based on PON deviation runs.

    For each gene, bins are scanned in genomic order and converted to a
    per-bin z-like deviation score:

        z_raw = (log2 - pon_log2) / pon_spread

    The z-scores are median-smoothed, then consecutive bins are grouped
    into runs when:
      - abs(z_smooth) >= config.z_bin_thresh
      - the sign of z_smooth stays constant within the run

    A run is promoted to a provisional gene region only if:
      - the gene has at least config.min_gene_targets bins
      - PON spread is available for the gene
      - the run contains at least config.min_run_bins bins
      - the run-level aggregated z-score is at least config.z_run_thresh

    Bins not assigned to any run keep the default region_id "no_region".

    Parameters
    ----------
    bins
        Per-bin table with at least:
            chr, gene.symbol, start, log2, pon_log2, pon_spread
    config
        Thresholds controlling run detection and smoothing.

    Returns
    -------
    pd.DataFrame
        Copy of `bins` with a provisional `region_id` column added/updated.
    """
    out = bins.copy()
    out["region_id"] = "no_region"

    next_region_id = 0

    def _run_abs_z(df_run: pd.DataFrame) -> float:
        """
        Compute absolute run-level z-score using mean effect and mean spread.
        """
        effect = df_run["log2"].to_numpy() - df_run["pon_log2"].fillna(0.0).to_numpy()
        if effect.size == 0:
            return np.nan

        mean_effect = float(np.nanmean(effect))
        if not np.isfinite(mean_effect):
            return np.nan

        spread = df_run["pon_spread"].fillna(0.0).to_numpy()
        spread_safe = np.where(spread <= 0, 1e-3, spread)
        mean_spread = float(np.nanmean(spread_safe)) if spread_safe.size else np.nan
        if not np.isfinite(mean_spread) or mean_spread <= 0:
            return np.nan

        sigma_effect = mean_spread / np.sqrt(float(len(effect)))
        if sigma_effect <= 0:
            return np.nan

        return abs(mean_effect) / sigma_effect

    def _label_gene_runs(gene_bins: pd.DataFrame) -> None:
        nonlocal next_region_id

        if gene_bins.shape[0] < config.min_gene_targets:
            return
        if gene_bins["pon_spread"].isna().all():
            return

        gene_bins = gene_bins.sort_values("start", kind="stable")

        bin_indices = gene_bins.index.to_numpy()

        effect = (
            gene_bins["log2"].to_numpy() - gene_bins["pon_log2"].fillna(0.0).to_numpy()
        )
        spread = gene_bins["pon_spread"].fillna(0.0).to_numpy()
        spread_safe = np.where(spread <= 0, 1e-3, spread)

        z_raw = effect / spread_safe

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

        current_run: list[int] = []
        current_sign: int | None = None

        def _finalize_run(run_indices: list[int]) -> None:
            nonlocal next_region_id

            if len(run_indices) < config.min_run_bins:
                return

            df_run = out.loc[run_indices]
            run_z = _run_abs_z(df_run)
            if not np.isfinite(run_z) or run_z < config.z_run_thresh:
                return

            region_label = f"generegion_{next_region_id}"
            next_region_id += 1
            out.loc[run_indices, "region_id"] = region_label

        for bin_index, z_value in zip(bin_indices, z_smooth):
            if not np.isfinite(z_value) or abs(z_value) < config.z_bin_thresh:
                if current_run:
                    _finalize_run(current_run)
                    current_run = []
                    current_sign = None
                continue

            z_sign = 1 if z_value > 0 else -1

            if current_sign is None or z_sign == current_sign:
                current_run.append(bin_index)
                current_sign = z_sign
            else:
                _finalize_run(current_run)
                current_run = [bin_index]
                current_sign = z_sign

        if current_run:
            _finalize_run(current_run)

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
         - the outer runs have similar mean effect
           (difference <= `config.bridge_delta`)

    2. Small-run merge:
       Merge a run with the previous run when:
         - their mean effects are similar
           (difference <= `config.merge_delta`)
         - at least one of the two runs has at most `config.small_seg_n` bins

    Effect is defined as:
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
              - mean_eff: mean(log2 - pon_log2) over the run
        eff_all
            Per-bin effect array aligned to gene_bins row order.
        """
        effect_all = (
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
                mean_effect = (
                    float(np.nanmean(effect_all[current_positions]))
                    if current_positions
                    else np.nan
                )
                runs.append(
                    {
                        "positions": current_positions,
                        "indices": run_indices,
                        "mean_eff": mean_effect,
                    }
                )
                current_label = region_labels[pos]
                current_positions = [pos]

        if current_positions:
            run_indices = gene_bins.index[current_positions].tolist()
            mean_effect = (
                float(np.nanmean(effect_all[current_positions]))
                if current_positions
                else np.nan
            )
            runs.append(
                {
                    "positions": current_positions,
                    "indices": run_indices,
                    "mean_eff": mean_effect,
                }
            )

        return runs, effect_all

    def _merge_gene_runs(gene_bins: pd.DataFrame) -> None:
        """
        Apply bridge-merge and small-run merge rules to one gene.
        """
        if gene_bins.empty or gene_bins.shape[0] <= 1:
            return

        gene_bins = gene_bins.sort_values("start", kind="stable")
        runs, effect_all = _collect_gene_runs(gene_bins)

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
                    and np.isfinite(run_a["mean_eff"])
                    and np.isfinite(run_c["mean_eff"])
                    and abs(run_a["mean_eff"] - run_c["mean_eff"])
                    <= config.bridge_delta
                ):
                    merged_positions = (
                        run_a["positions"] + run_b["positions"] + run_c["positions"]
                    )
                    merged_indices = (
                        run_a["indices"] + run_b["indices"] + run_c["indices"]
                    )
                    merged_mean_eff = (
                        float(np.nanmean(effect_all[merged_positions]))
                        if merged_positions
                        else np.nan
                    )

                    merged_runs.append(
                        {
                            "positions": merged_positions,
                            "indices": merged_indices,
                            "mean_eff": merged_mean_eff,
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
                    np.isfinite(current_run["mean_eff"])
                    and np.isfinite(previous_run["mean_eff"])
                    and abs(current_run["mean_eff"] - previous_run["mean_eff"])
                    <= config.merge_delta
                    and (
                        len(current_run["indices"]) <= config.small_seg_n
                        or len(previous_run["indices"]) <= config.small_seg_n
                    )
                )

                if should_merge:
                    merged_positions = (
                        previous_run["positions"] + current_run["positions"]
                    )
                    merged_indices = previous_run["indices"] + current_run["indices"]
                    merged_mean_eff = (
                        float(np.nanmean(effect_all[merged_positions]))
                        if merged_positions
                        else np.nan
                    )

                    merged_runs[-1] = {
                        "positions": merged_positions,
                        "indices": merged_indices,
                        "mean_eff": merged_mean_eff,
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

    Parameters
    ----------
    bins
        Per-bin table containing at least:
            chr, gene.symbol, region_id, start, end, log2, pon_log2, pon_spread

    Returns
    -------
    pd.DataFrame
        One row per (chr, gene.symbol, region_id) with columns including:
            region_start, region_end, n_targets, mean_log2,
            min_log2, max_log2, pon_mean_log2, pon_mean_spread, n.targets
    """
    agg_dict: dict[str, list[str]] = {
        "start": ["min"],
        "end": ["max"],
        "log2": ["count", "mean", "min", "max"],
        "pon_log2": ["mean"],
        "pon_spread": ["mean"],
    }

    regions_df = (
        bins.groupby(["chr", "gene.symbol", "region_id"], as_index=False)
        .agg(agg_dict)
        .copy()
    )

    regions_df = _flatten_agg_columns(regions_df).rename(
        columns={
            "start_min": "region_start",
            "end_max": "region_end",
            "log2_count": "n_targets",
            "log2_mean": "mean_log2",
            "log2_min": "min_log2",
            "log2_max": "max_log2",
            "pon_log2_mean": "pon_mean_log2",
            "pon_spread_mean": "pon_mean_spread",
        }
    )

    # Keep compatibility with older downstream code using n.targets
    regions_df["n.targets"] = regions_df["n_targets"]

    return regions_df


def _score_pon_regions(
    regions_df: pd.DataFrame,
    *,
    config: GeneRegionConfig,
) -> pd.DataFrame:
    """
    Compute PON-based effect, z-score, signal class, and CNV indication
    for each gene-region.

    Scoring is based on deviation from the PON baseline:

        pon_region_effect = mean_log2 - pon_mean_log2

    The absolute effect is converted into a z-like score using the mean
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
            pon_region_effect
            pon_region_z
            pon_region_direction
            pon_region_signal
            pon_region_indication
    """
    out = regions_df.copy()

    # Mean deviation from the PON baseline
    out["pon_region_effect"] = out["mean_log2"] - out["pon_mean_log2"]

    # Region-level z-score using target count as effective sample size
    out["pon_region_z"] = out.apply(
        lambda row: _pon_abs_z(
            row["pon_region_effect"],
            row["pon_mean_spread"],
            row.get("n_targets", row.get("n.targets", np.nan)),
            min_n=config.min_region_targets_for_call,
        ),
        axis=1,
    )

    # Direction of deviation relative to PON baseline
    out["pon_region_direction"] = out["pon_region_effect"].apply(_pon_direction)

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
        lambda row: _pon_cnv_call_from_effect(
            is_strong=str(row.get("pon_region_signal", "")).strip().lower() == "strong",
            effect_log2=row.get("pon_region_effect", np.nan),
            gain_gt=config.pon_gain_gt,
            loss_lt=config.pon_loss_lt,
            weak_value="NEUTRAL",
            neutral_value="NEUTRAL",
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
    regions_df = add_cnv_calls_wide(regions_df, sex=sex)

    # Drop columns if present
    drop_cols = [
        c
        for c in [
            "region_id",
            "n_targets",
            "pon_region_direction",
            "cnvkit_seg_start",
            "cnvkit_seg_end",
            "cnvkit_seg_log2",
            "purecn_seg_start",
            "purecn_seg_end",
            "purecn_seg_mean_log2",
            "purecn_num_snps",
            "purecn_maf_observed",
            "cnvkit_seg_cn1",
            "cnvkit_seg_cn2",
            "purecn_M",
            "purecn_M_flagged",
            "pon_region_effect",
            "cnvkit_seg_depth",
        ]
        if c in regions_df.columns
    ]
    if drop_cols:
        regions_df = regions_df.drop(columns=drop_cols)

    return reorder_and_sort_table(regions_df, GENE_TABLE_SPEC)


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
    segments = add_cnv_calls_wide(
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
    return reorder_and_sort_table(segments, SEGMENT_TABLE_SPEC)
