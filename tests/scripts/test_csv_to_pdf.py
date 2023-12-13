"""Test converting CSV to PDF."""
from pathlib import Path

from click.testing import CliRunner, Result
from pypdf import PdfReader

from BALSAMIC.assets.scripts.csv_to_pdf import csv_to_pdf, get_table_html_page


def test_csv_to_pdf(
    purity_csv_path: Path, tmp_path: Path, cli_runner: CliRunner
) -> None:
    """Test converting of a CSV file to PDF."""

    # GIVEN an input CSV file

    # GIVEN an output PDF file
    pdf_path: Path = Path(tmp_path, "csv_to_pdf.pdf")

    # WHEN converting the CSV file to PDF
    result: Result = cli_runner.invoke(
        csv_to_pdf, [purity_csv_path.as_posix(), pdf_path.as_posix()]
    )

    # THEN the output PDF file should exist
    assert result.exit_code == 0
    assert pdf_path.is_file()

    # THEN the output PDF file should contain the CSV table
    reader: PdfReader = PdfReader(stream=pdf_path)
    pdf_page: str = reader.pages[0].extract_text()
    assert purity_csv_path.stem in pdf_page


def test_get_table_html_page():
    """Test table insertion in an HTML page"""

    # GIVEN an HTML table and a test table name
    html_table: str = """
        <table>
          <tr>
            <th>Header 1</th>
            <th>Header 2</th>
          </tr>
          <tr>
            <td>Data 1</td>
            <td>Data 2</td>
          </tr>
        </table>
    """
    table_name: str = "Test Table"

    # WHEN adding the table to an HTML page
    html_page: str = get_table_html_page(html_table=html_table, table_name=table_name)

    # THEN the table HTML page should be successfully created
    assert "<html>" in html_page
    assert f"<h2>{table_name}</h2>" in html_page
    assert html_table in html_page
