from __future__ import annotations

from pathlib import Path

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
        "Curated",
        "Comment",
    ]
    keep = [c for c in wanted_cols if c in df.columns]
    return df[keep] if keep else df


def compute_dlr_spread_from_cnr(
    cnr_df: pd.DataFrame,
    exclude_sex_chromosomes: bool = True,
) -> float:
    """
    Compute a DLR-like spread metric from a CNVkit .cnr file.

    - Sort bins by chr, start.
    - Compute differences between adjacent log2 values.
    - Return np.nanstd of these differences.

    Returns np.nan if fewer than 3 usable bins.
    """

    cnr_df["type"] = np.where(cnr_df["gene.symbol"] == "Antitarget", "Antitarget", "Target")

    cnr_df = cnr_df[cnr_df["type"] == "Target"]

    if exclude_sex_chromosomes:
        cnr_df = cnr_df[~cnr_df["chr"].isin(["X", "Y"])]

    cnr_df = cnr_df.dropna(subset=["log2"])
    if cnr_df.empty:
        return float("nan")

    cnr_df["chr_sort"] = cnr_df["chr"].map(chrom_sort_key)
    cnr_df = cnr_df.sort_values(["chr_sort", "start"]).reset_index(drop=True)
    log2_vals = cnr_df["log2"].to_numpy()

    if log2_vals.size < 3:
        return float("nan")

    diffs = np.diff(log2_vals)
    dlr = float(np.nanstd(diffs))
    return dlr


def compute_pon_spread_summaries(
    pon_df: pd.DataFrame,
    exclude_sex_chromosomes: bool = True,
) -> dict[str, float]:
    """
    Compute simple spread summaries from a CNVkit PON .cnn file.

    Returns a dict with:
      - pon_spread_median_target
      - pon_spread_q90_target

    (np.nan if not available)
    """

    chr_col = "chr"

    # mark Target / Antitarget
    pon_df["type"] = np.where(pon_df["gene.symbol"] == "Antitarget", "Antitarget", "Target")

    pon_df = pon_df[pon_df["type"] == "Target"]

    if exclude_sex_chromosomes:
        pon_df = pon_df[~pon_df[chr_col].isin(["X", "Y"])]

    pon_df = pon_df.dropna(subset=["pon_spread"])
    if pon_df.empty:
        return {
            "pon_spread_median_target": float("nan"),
            "pon_spread_q90_target": float("nan"),
        }

    return {
        "pon_spread_median_target": float(pon_df["pon_spread"].median()),
        "pon_spread_q90_target": float(pon_df["pon_spread"].quantile(0.90)),
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
    metrics["Log2-spread"] = compute_dlr_spread_from_cnr(cnr_df)

    # 2) PON spread summaries
    if pon_df is not None:
        pon_stats = compute_pon_spread_summaries(pon_df)
        metrics.update(pon_stats)
    else:
        metrics["pon_spread_median_target"] = float("nan")
        metrics["pon_spread_q90_target"] = float("nan")

    return pd.DataFrame([metrics])
