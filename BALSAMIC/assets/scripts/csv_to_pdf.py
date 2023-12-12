"""Script for converting a CSV file to a PDF."""
from pathlib import Path

import click
from pandas import DataFrame, read_csv

from BALSAMIC.utils.pdf_report import get_table_html, html_to_pdf


@click.command()
@click.argument(
    "csv_path", nargs=1, required=True, type=click.Path(exists=True, resolve_path=True)
)
@click.argument("pdf_path", nargs=1, required=True, type=click.Path(resolve_path=True))
@click.option(
    "--delimiter",
    type=click.STRING,
    default=",",
    show_default=True,
    help="CSV file delimiter",
)
def csv_to_pdf(csv_path: str, pdf_path: str, delimiter: str) -> None:
    """Convert CSV file to a PDF."""
    df: DataFrame = read_csv(filepath_or_buffer=csv_path, delimiter=delimiter)
    html_table: str = df.to_html(
        index=False, na_rep="NA", justify="center", escape=False
    )
    html_page: str = get_table_html(
        html_table=html_table, table_name=Path(csv_path).stem
    )
    html_to_pdf(html_string=html_page, pdf_path=pdf_path)


if __name__ == "__main__":
    csv_to_pdf()
