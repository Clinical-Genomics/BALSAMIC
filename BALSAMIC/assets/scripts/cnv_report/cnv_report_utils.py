from __future__ import annotations

# Third-party
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
