from __future__ import annotations

from pathlib import Path

import math
import numpy as np
import pandas as pd


from cnv_report_utils import chrom_sort_key


def read_purecn_summary(purity_csv: str | Path | None) -> pd.DataFrame | None:
    """Read and return PureCN purity table"""
    df = pd.read_csv(purity_csv)
    wanted_cols = [
        "Sampleid",
        "Purity",
        "Ploidy",
        "Sex",
        "Contamination",
        "Flagged",
        "Failed",
        "Comment",
    ]
    keep = [c for c in wanted_cols if c in df.columns]
    return df[keep] if keep else df


def compute_dlr_spread_from_cnr(
    cnr_df: pd.DataFrame,
    exclude_sex_chromosomes: bool = True,
) -> float:
    """
    Compute a DLRSpread-like CNV noise metric from a CNVkit .cnr-style dataframe.

    Steps:
      - Keep Target bins only
      - Optionally exclude X/Y
      - Sort bins by chromosome and position
      - Compute log2 differences between adjacent bins *within each chromosome*
      - Return std(diff) / sqrt(2)

    Returns np.nan if too few usable bins/differences are available.
    """
    df = cnr_df.copy()

    df = df[df["gene.symbol"] != "Antitarget"].copy()

    if exclude_sex_chromosomes:
        df = df[~df["chr"].isin(["X", "Y"])].copy()

    df = df.dropna(subset=["chr", "start", "log2"])
    if df.empty:
        return float("nan")

    df["chr_sort"] = df["chr"].map(chrom_sort_key)
    df = df.sort_values(["chr_sort", "start"], kind="stable")

    diffs = df.groupby("chr", sort=False)["log2"].diff().dropna().to_numpy()

    if diffs.size < 2:
        return float("nan")

    return float(np.nanstd(diffs) / math.sqrt(2))


def compute_filtered_out_bin_metrics(
    cnr_df: pd.DataFrame,
    pon_df: pd.DataFrame,
    *,
    key_cols: tuple[str, str, str] = ("chr", "start", "end"),
) -> dict[str, float | int]:
    """
    Summarize bins present in PON but missing from CNR ("filtered out").

    All calculations use UNIQUE (chr, start, end) bins.
    """

    # Unique bin definitions
    cnr_keys = cnr_df[list(key_cols)].drop_duplicates()
    pon_unique = pon_df.drop_duplicates(subset=list(key_cols)).copy()
    pon_keys = pon_unique[list(key_cols)]

    # Missing (filtered out) bins
    missing_keys = (
        pon_keys.merge(cnr_keys, on=list(key_cols), how="left", indicator=True)
        .query('_merge == "left_only"')
        .drop(columns="_merge")
    )

    n_total = int(len(pon_unique))
    n_filtered = int(len(missing_keys))

    pct_filtered = (
        round(float(n_filtered / n_total * 100.0), 3) if n_total else float("nan")
    )

    return {
        "Targets total": n_total,
        "Targets filtered": f"{n_filtered} ({pct_filtered}%)",
    }


def compute_summary_metrics(
    cnr_df: pd.DataFrame,
    pon_df: pd.DataFrame | None,
) -> pd.DataFrame:
    """
    Build a 1-row DataFrame with QC / burden summaries.

    Includes:
      - Log2-spread (DLR-like) from the CNVkit .cnr file
      - PON spread summaries from the CNVkit .cnn file (targets only)

    Parameters
    ----------
    cnr_df
    pon_df


    Returns
    -------
    pd.DataFrame
        Single-row DataFrame of metrics.
    """
    metrics: dict[str, float | int] = {}

    # 1) DLR-like spread
    metrics["Log2-noise"] = compute_dlr_spread_from_cnr(cnr_df)

    # 2) Add dropped bins summary
    if pon_df is not None:
        metrics.update(compute_filtered_out_bin_metrics(cnr_df, pon_df))

    return pd.DataFrame([metrics])
