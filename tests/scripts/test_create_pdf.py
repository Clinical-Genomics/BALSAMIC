from pathlib import Path

from BALSAMIC.assets.scripts.generate_cnv_report import (
    get_pdf_instance,
    add_data_to_pdf,
    add_plots_to_pdf,
    generate_cnv_report,
    PDF,
    get_pdf_data,
)


def test_get_pdf_data():
    """Test pdf data extraction from a list of files."""

    # GIVEN a list of input files
    data_paths = ["statistics.txt", "plot_0.png", "plot_1.png"]

    # WHEN retrieving the statistics and plots tuple
    statistics, plots = get_pdf_data(data_paths)

    # THEN the expected files should be returned
    assert data_paths[0] in statistics
    assert data_paths[1] in plots
    assert data_paths[2] in plots


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
    statistics_paths = [
        "tests/test_data/cnv_report/CNV.somatic.sample_tumor_normal_wgs.ascat.samplestatistics.txt"
    ]

    # WHEN generating the PDF with the statistics
    pdf: PDF = add_data_to_pdf(pdf=pdf, data_paths=statistics_paths)

    # THEN check if the statistics are appended to the created PDF
    assert isinstance(pdf, PDF)
    assert pdf.page_no() == 1


def test_add_plots_to_pdf():
    """Test plots appending to a PDF file."""

    # GIVEN a PDF instance and some dummy PNG plots
    pdf: PDF = get_pdf_instance()
    plot_paths = [
        "tests/test_data/cnv_report/CNV.somatic.sample_tumor_normal_wgs.ascat.sunrise.png",
        "tests/test_data/cnv_report/CNV.somatic.sample_tumor_normal_wgs.ascat.germline.png",
    ]

    # WHEN adding the plots to a PDF instance
    pdf: PDF = add_plots_to_pdf(pdf, plot_paths)

    # THEN check if the images are correctly appended to the PDF
    assert isinstance(pdf, PDF)
    assert pdf.page_no() == len(plot_paths)


def test_generate_cnv_report_tumor_normal(tmp_path, cli_runner):
    """Test generation of a PDF report for a WGS TN case."""

    # GIVEN dummy input data and plots
    data_paths = [
        "tests/test_data/cnv_report/CNV.somatic.sample_tumor_normal_wgs.ascat.samplestatistics.txt",
        "tests/test_data/cnv_report/CNV.somatic.sample_tumor_normal_wgs.ascat.germline.png",
        "tests/test_data/cnv_report/CNV.somatic.sample_tumor_normal_wgs.ascat.sunrise.png",
        "tests/test_data/cnv_report/CNV.somatic.sample_tumor_only_wgs.cnvpytor.circular.png",
        "tests/test_data/cnv_report/CNV.somatic.sample_tumor_only_wgs.cnvpytor.scatter.png",
    ]

    # GIVEN the output path
    output_path: Path = Path(tmp_path, "report.pdf")

    # WHEN invoking the python script
    result = cli_runner.invoke(
        generate_cnv_report,
        [
            data_paths[0],
            data_paths[1],
            data_paths[2],
            data_paths[3],
            data_paths[4],
            "--output",
            output_path,
        ],
    )

    # THEN check if the PDF is correctly created and there is no errors
    assert result.exit_code == 0
    assert Path(output_path).exists()


def test_generate_cnv_report_tumor_only(tmp_path, cli_runner):
    """Test generation of a PDF report for a WGS TO case."""

    # GIVEN dummy input data and plots
    plot_paths = [
        "tests/test_data/cnv_report/CNV.somatic.sample_tumor_only_wgs.cnvpytor.circular.png",
        "tests/test_data/cnv_report/CNV.somatic.sample_tumor_only_wgs.cnvpytor.scatter.png",
    ]

    # GIVEN the output path
    output_path: Path = Path(tmp_path, "report.pdf")

    # WHEN invoking the python script
    result = cli_runner.invoke(
        generate_cnv_report,
        [
            plot_paths[0],
            plot_paths[1],
            "--output",
            output_path,
        ],
    )

    # THEN check if the PDF is correctly created and there is no errors
    assert result.exit_code == 0
    assert Path(output_path).exists()
