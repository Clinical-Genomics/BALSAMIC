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
    cnr_path: str | Path,
    exclude_sex_chromosomes: bool = True,
) -> float:
    """
    Compute a DLR-like spread metric from a CNVkit .cnr file.

    - Sort bins by chr, start.
    - Compute differences between adjacent log2 values.
    - Return np.nanstd of these differences.

    Returns np.nan if fewer than 3 usable bins.
    """
    cnr = pd.read_csv(cnr_path, sep="\t")

    chr_col = "chromosome"

    cnr["type"] = np.where(cnr["gene"] == "Antitarget", "Antitarget", "Target")

    cnr = cnr[cnr["type"] == "Target"]

    if exclude_sex_chromosomes:
        cnr = cnr[~cnr[chr_col].isin(["X", "Y"])]

    cnr = cnr.dropna(subset=["log2"])
    if cnr.empty:
        return float("nan")

    cnr["chr_sort"] = cnr[chr_col].map(chrom_sort_key)
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

    chr_col = "chromosome"

    # mark Target / Antitarget
    cnn["type"] = np.where(cnn["gene"] == "Antitarget", "Antitarget", "Target")

    cnn = cnn[cnn["type"] == "Target"]

    if exclude_sex_chromosomes:
        cnn = cnn[~cnn[chr_col].isin(["X", "Y"])]

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

    return pd.DataFrame([metrics])
