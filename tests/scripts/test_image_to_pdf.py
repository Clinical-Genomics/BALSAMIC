"""Test converting images to PDF."""
from pathlib import Path

from click.testing import CliRunner, Result
from pypdf import PdfReader

from BALSAMIC.assets.scripts.image_to_pdf import image_to_pdf


def test_image_to_pdf(
    cnv_plot_path: Path, tmp_path: Path, cli_runner: CliRunner
) -> None:
    """Test converting of an image file to PDF."""

    # GIVEN an input CNV plot file

    # GIVEN an output PDF file
    pdf_path: Path = Path(tmp_path, "image_to_pdf.pdf")

    # WHEN converting the plot to PDF
    result: Result = cli_runner.invoke(
        image_to_pdf, [cnv_plot_path.as_posix(), pdf_path.as_posix()]
    )

    # THEN the output PDF file should exist
    assert result.exit_code == 0
    assert pdf_path.is_file()

    # THEN the output PDF file should contain the image
    reader: PdfReader = PdfReader(stream=pdf_path)
    pdf_page: str = reader.pages[0].extract_text()
    assert cnv_plot_path.stem in pdf_page
