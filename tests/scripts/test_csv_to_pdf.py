"""Test converting CSV to PDF."""
from pathlib import Path

from click.testing import CliRunner, Result
from pypdf import PdfReader

from BALSAMIC.assets.scripts.csv_to_pdf import csv_to_pdf


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
