"""Test PDF file merging."""
from pathlib import Path
from typing import List

from click.testing import CliRunner, Result
from pypdf import PdfReader, PdfWriter

from BALSAMIC.assets.scripts.merge_pdfs import merge_pdfs


def test_merge_pdfs(tmp_path: Path, cli_runner: CliRunner) -> None:
    """Test merging of multiple PDF files."""

    # GIVEN a list of empty PDF files
    input_pdfs: List[str] = [
        Path(tmp_path, f"file_{pdf}.pdf").as_posix() for pdf in [1, 2, 3]
    ]
    for pdf in input_pdfs:
        pdf_writer = PdfWriter()
        pdf_writer.add_blank_page(111, 111)
        pdf_writer.write(pdf)

    # GIVEN an output PDF file and PDF reader
    output_pdf: str = Path(tmp_path, "output.pdf").as_posix()

    # WHEN merging multiple PDFs
    result: Result = cli_runner.invoke(merge_pdfs, input_pdfs + [output_pdf])

    # THEN the output file should contain all the PDF files
    pdf_reader = PdfReader(output_pdf)
    assert len(pdf_reader.pages) == len(input_pdfs)
    assert result.exit_code == 0
