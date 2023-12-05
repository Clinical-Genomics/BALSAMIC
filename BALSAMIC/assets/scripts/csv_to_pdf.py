"""Script for converting a CSV file to a PDF."""
from pathlib import Path

import click
import pdfkit
from pandas import DataFrame, read_csv


@click.command()
@click.argument(
    "csv_path", nargs=1, required=True, type=click.Path(exists=True, resolve_path=True)
)
@click.argument("pdf_path", nargs=1, required=True, type=click.Path(resolve_path=True))
def csv_to_pdf(csv_path: str, pdf_path: str) -> None:
    """Convert CSV file to a PDF."""
    df: DataFrame = read_csv(csv_path)
    html_table: str = df.to_html(
        index=False, na_rep="NA", justify="center", escape=False
    )
    html_page: str = get_table_html_page(
        html_table=html_table, table_name=Path(csv_path).stem
    )
    pdfkit.from_string(
        input=html_page,
        output_path=pdf_path,
        options={
            "page-size": "A4",
            "orientation": "landscape",
            "enable-local-file-access": None,
        },
    )


def get_table_html_page(html_table: str, table_name: str) -> str:
    """Return HTML-rendered content with the provided HTML table."""
    return f"""
        <!DOCTYPE html>
        <html>
        <head>
            <style>
                h2 {{text-align: center; padding: 10px;}}
                table {{margin: 0 auto; border: 1px solid black; border-collapse: collapse; text-align: center;}}
                th {{font-size: 12pt; padding: 5px; background: #cccccc;}}
                td {{font-size: 10pt; padding: 5px;}}
                tr:nth-child(even) {{background: #eeeeee;}}
            </style>
        </head>
        <body>
            <h2>{table_name}</h2>
            {html_table}
        </body>
        </html>
    """


if __name__ == "__main__":
    csv_to_pdf()
