"""Script to merge PDFs."""
from typing import List

import click
from pypdf import PdfWriter


@click.command()
@click.argument(
    "input_pdfs",
    nargs=-1,
    required=True,
    type=click.Path(exists=True, resolve_path=True),
)
@click.argument(
    "output_pdf", nargs=1, required=True, type=click.Path(resolve_path=True)
)
def merge_pdfs(input_pdfs: List[str], output_pdf: str):
    """Merge PDFs into a single file."""
    merger: PdfWriter = PdfWriter()
    for pdf in input_pdfs:
        merger.append(pdf)
    merger.write(output_pdf)
    merger.close()


if __name__ == "__main__":
    merge_pdfs()
