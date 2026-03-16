from __future__ import annotations

# Standard library

from pathlib import Path

# Third-party
import pandas as pd
import fitz


def pdf_first_page_to_png(pdf_path: str, png_path: str, dpi: int = 300) -> None:
    """
    Render the first page of a PDF to a PNG image using PyMuPDF.
    """
    doc = fitz.open(str(pdf_path))
    try:
        page = doc[0]
        page.get_pixmap(dpi=dpi).save(str(png_path))
    finally:
        doc.close()


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



def _left_merge_pon(cnr_bins: pd.DataFrame, pon_bins: pd.DataFrame) -> pd.DataFrame:
    """Left-merge PON columns onto bins by (chr,start,end)."""
    if pon_bins is None or pon_bins.empty:
        return cnr_bins
    return cnr_bins.merge(pon_bins, how="left", on=["chr", "start", "end"])
