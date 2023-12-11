"""Test utility function for PDF generation."""
from pathlib import Path

from pypdf import PdfReader

from BALSAMIC.utils.pdf_report import get_table_html_page, html_to_pdf


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


def test_html_to_pdf(tmp_path: Path):
    """Test PDF file generation from HTML string."""

    # GIVEN an HTML string
    html_string: str = """
        <html>
        <body>
            <p>Hello!</p>
        </body>
        </html>
    """

    # GIVEN an output PDF file
    pdf_path: Path = Path(tmp_path, "csv_to_pdf.pdf")

    # WHEN generating the pdf file
    html_to_pdf(html_string=html_string, pdf_path=pdf_path)

    # THEN the output PDF file should exist
    assert pdf_path.is_file()

    # THEN the output PDF file should contain the mock HTML string
    reader: PdfReader = PdfReader(stream=pdf_path)
    pdf_page: str = reader.pages[0].extract_text()
    assert "Hello!" in pdf_page
