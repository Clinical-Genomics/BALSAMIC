from __future__ import annotations

# Standard library
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Any, Mapping, Sequence

# Third-party
import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError
import fitz

# Local
from BALSAMIC.constants.analysis import Gender


# =============================================================================
# Generic helpers
# =============================================================================


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


def strip_chr_prefix(series: pd.Series) -> pd.Series:
    """Normalize chromosome values by stripping a leading 'chr' prefix."""
    return series.astype(str).str.replace("^chr", "", regex=True)


def chrom_sort_key(chrom: str) -> tuple[int, int | str]:
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


def stable_sort_by_chr_interval(
    df: pd.DataFrame,
    chr_col: str,
    start_col: str,
    end_col: str,
    *,
    tmp_col: str = "chr_sort",
) -> pd.DataFrame:
    """Stable-sort by (chr, start, end) using `_chrom_sort_key`."""
    df = df.copy()
    df[tmp_col] = df[chr_col].map(chrom_sort_key)
    df = df.sort_values(by=[tmp_col, start_col, end_col], kind="stable").drop(
        columns=[tmp_col]
    )
    return df


def reorder_columns(df: pd.DataFrame, preferred: list[str]) -> pd.DataFrame:
    """Move preferred columns first (when present), preserving the rest."""
    preferred_present = [c for c in preferred if c in df.columns]
    return df[preferred_present + [c for c in df.columns if c not in preferred_present]]


def detect_chr_col(df: pd.DataFrame, candidates: list[str]) -> str:
    for c in candidates:
        if c in df.columns:
            return c
    raise ValueError(f"Could not find chromosome column among: {candidates}")


def flatten_agg_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Flatten MultiIndex columns produced by pandas .agg()."""
    df = df.copy()
    df.columns = [
        "_".join(col).strip("_") if isinstance(col, tuple) else col
        for col in df.columns
    ]
    return df


def explode_multigene_bins(cnr: pd.DataFrame) -> pd.DataFrame:
    """Drop pseudo genes and explode multi-gene bins into one row per gene.symbol."""
    cnr = cnr.copy()

    gene_col = "gene"

    # drop NA and pseudo-genes
    cnr = cnr[cnr[gene_col].notna()]
    cnr = cnr[~cnr[gene_col].isin(["Antitarget", "-"])]

    # split on commas, strip whitespace
    gene_series = cnr[gene_col].astype(str).str.split(r"\s*,\s*", regex=True)

    cnr = cnr.assign(gene_symbol=gene_series).explode("gene_symbol")
    cnr = cnr[
        cnr["gene_symbol"].notna() & (cnr["gene_symbol"].astype(str).str.len() > 0)
    ]
    cnr = cnr.rename(columns={"gene_symbol": "gene.symbol"})
    return cnr


def _left_merge_pon(cnr_bins: pd.DataFrame, pon_bins: pd.DataFrame) -> pd.DataFrame:
    """Left-merge PON columns onto bins by (chr,start,end)."""
    if pon_bins is None or pon_bins.empty:
        return cnr_bins
    return cnr_bins.merge(pon_bins, how="left", on=["chr", "start", "end"])
