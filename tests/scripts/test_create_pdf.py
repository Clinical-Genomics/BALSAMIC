from pathlib import Path

from BALSAMIC.assets.scripts.generate_cnv_report import (
    get_pdf_instance,
    add_data_to_pdf,
    add_plots_to_pdf,
    generate_cnv_report,
    PDF,
)


def test_get_pdf_instance():
    """Test FPDF instance generation."""

    # WHEN creating a dummy FPDF file
    pdf: PDF = get_pdf_instance()

    # THEN check if the PDF has been correctly created
    assert isinstance(pdf, PDF)


def test_add_data_to_pdf():
    """Test add statistics to a PDF instance."""

    # GIVEN a PDF instance and an output sample statistics .txt file
    pdf: PDF = get_pdf_instance()
    statistics_path = (
        "tests/test_data/ascat_output/CNV.somatic.SAMPLE.ascat.samplestatistics.txt"
    )

    # WHEN generating the PDF with the statistics
    pdf: PDF = add_data_to_pdf(pdf=pdf, data_path=statistics_path)

    # THEN check if the statistics are appended to the created PDF
    assert isinstance(pdf, PDF)
    assert pdf.page_no() == 1


def test_add_plots_to_pdf():
    """Test plots appending to a PDF file."""

    # GIVEN a PDF instance and some dummy PNG plots
    pdf: PDF = get_pdf_instance()
    plot_paths = [
        "tests/test_data/ascat_output/CNV.somatic.SAMPLE.ascat.sunrise.png",
        "tests/test_data/ascat_output/CNV.somatic.SAMPLE.ascat.germline.png",
    ]

    # WHEN adding the plots to a PDF instance
    pdf: PDF = add_plots_to_pdf(pdf, plot_paths)

    # THEN check if the images are correctly appended to the PDF
    assert isinstance(pdf, PDF)
    assert pdf.page_no() == len(plot_paths)


def test_generate_cnv_report(tmp_path, cli_runner):
    """Test generation of a PDF report."""

    # GIVEN dummy input data and plots
    statistics_path = (
        "tests/test_data/ascat_output/CNV.somatic.SAMPLE.ascat.samplestatistics.txt"
    )
    plot_paths = [
        "tests/test_data/ascat_output/CNV.somatic.SAMPLE.ascat.germline.png",
        "tests/test_data/ascat_output/CNV.somatic.SAMPLE.ascat.sunrise.png",
    ]

    # GIVEN the output path
    output_path: Path = Path(tmp_path, "ascat.output.pdf")

    # WHEN invoking the python script
    result = cli_runner.invoke(
        generate_cnv_report,
        [
            statistics_path,
            plot_paths[0],
            plot_paths[1],
            "--output",
            output_path,
        ],
    )

    # THEN check if the PDF is correctly created and there is no errors
    assert result.exit_code == 0
    assert Path(output_path).exists()
