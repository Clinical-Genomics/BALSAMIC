"""Script for converting images to PDF."""
from pathlib import Path

import click

from BALSAMIC.utils.pdf_report import get_image_html, html_to_pdf


@click.command()
@click.argument(
    "image_path",
    nargs=1,
    required=True,
    type=click.Path(exists=True, resolve_path=True),
)
@click.argument("pdf_path", nargs=1, required=True, type=click.Path(resolve_path=True))
def image_to_pdf(image_path: str, pdf_path: str) -> None:
    """Convert image file to a PDF."""
    html_page: str = get_image_html(
        image_path=Path(image_path), image_name=Path(image_path).stem
    )
    html_to_pdf(html_string=html_page, pdf_path=pdf_path)


if __name__ == "__main__":
    image_to_pdf()
