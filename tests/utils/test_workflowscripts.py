from unittest import mock
from pathlib import Path
import pytest
from fpdf import FPDF

from BALSAMIC.utils.cli import generate_h5
from BALSAMIC.utils.workflowscripts import (
    plot_analysis,
    create_pdf,
    add_images_pdf,
    add_table_pdf,
    save_ascat_output_pdf,
)


def test_plot_analysis(tmp_path_factory):
    # GIVEN a dummy log file
    dummy_log_file = Path(
        "tests/test_data/dummy_run_logs/BALSAMIC.T_panel.bwa_mem.123.sh_31415926535.err"
    )
    dummy_h5 = "tests/test_data/dummy_run_logs/BALSAMIC.T_panel.bwa_mem.123.h5"
    dummy_path = tmp_path_factory.mktemp("dummy_pdf_path")
    dummy_pdf_name = dummy_path / "BALSAMIC.T_panel.bwa_mem.123.pdf"
    dummy_pdf_name.touch()

    # WHEN calling plot_analysis
    actual_pdf_file = plot_analysis(dummy_log_file, dummy_h5, dummy_pdf_name)

    assert Path(actual_pdf_file).exists()


def test_plot_analysis_bad_h5(tmp_path_factory):
    # GIVEN a dummy log file
    dummy_log_file = Path(
        "tests/test_data/dummy_run_logs/BALSAMIC.T_panel.bwa_mem.123.sh_31415926535.err"
    )
    dummy_h5 = "tests/test_data/dummy_run_logs/bad_format.h5"
    dummy_path = tmp_path_factory.mktemp("dummy_pdf_path")
    dummy_pdf_name = dummy_path / "plot_file.pdf"
    dummy_pdf_name.touch()

    # WHEN calling plot_analysis
    actual_pdf_file = plot_analysis(dummy_log_file, dummy_h5, dummy_pdf_name)

    assert actual_pdf_file is None


def test_create_pdf():
    # WHEN creating a dummy FPDF file
    pdf = create_pdf()

    # THEN check if the pdf has been correctly created
    assert isinstance(pdf, FPDF)


def test_add_images_pdf():
    # GIVEN ascatNGgs output PNG images
    test_images_path = [
        "tests/test_data/ascat_output/CNV.somatic.SAMPLE.ascat.sunrise.png",
        "tests/test_data/ascat_output/CNV.somatic.SAMPLE.ascat.germline.png",
    ]

    # WHEN calling the function
    pdf = add_images_pdf(create_pdf(), test_images_path)

    # THEN check if the images are appended to the PDF
    assert isinstance(pdf, FPDF)
    assert pdf.page_no() == 2


def test_add_table_pdf():
    # GIVEN ascatNGgs output sample statistics .txt
    test_statistics_path = (
        "tests/test_data/ascat_output/CNV.somatic.SAMPLE.ascat.samplestatistics.txt"
    )

    # WHEN calling the function
    pdf = add_table_pdf(create_pdf(), test_statistics_path)

    # THEN check if the table is appended to the created PDF
    assert isinstance(pdf, FPDF)
    assert pdf.page_no() == 1


def test_save_ascat_output_pdf(tmp_path):
    # GIVEN ascatNGgs output files
    test_paths = [
        "tests/test_data/ascat_output/CNV.somatic.SAMPLE.ascat.samplestatistics.txt",
        "tests/test_data/ascat_output/CNV.somatic.SAMPLE.ascat.sunrise.png",
        "tests/test_data/ascat_output/CNV.somatic.SAMPLE.ascat.germline.png",
    ]

    # GIVEN the output path
    output_path = tmp_path / "ascat.output.pdf"

    # WHEN calling the function
    save_ascat_output_pdf(output_path, test_paths[0], test_paths[1:])

    # THEN check if the PDF is correctly created
    assert Path(output_path).exists()
