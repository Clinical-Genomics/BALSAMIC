from pathlib import Path

from fpdf import FPDF

from BALSAMIC.assets.scripts.create_pdf import (
    generate_fpdf,
    add_images_pdf,
    add_table_pdf,
    create_pdf,
)


def test_generate_fpdf():
    # WHEN creating a dummy FPDF file
    pdf = generate_fpdf()

    # THEN check if the pdf has been correctly created
    assert isinstance(pdf, FPDF)


def test_add_images_pdf():
    # GIVEN ascatNGgs output PNG images
    test_images_path = [
        "tests/test_data/ascat_output/CNV.somatic.SAMPLE.ascat.sunrise.png",
        "tests/test_data/ascat_output/CNV.somatic.SAMPLE.ascat.germline.png",
    ]

    # WHEN calling the function
    pdf = add_images_pdf(generate_fpdf(), test_images_path)

    # THEN check if the images are appended to the PDF
    assert isinstance(pdf, FPDF)
    assert pdf.page_no() == 2


def test_add_table_pdf():
    # GIVEN ascatNGgs output sample statistics .txt
    test_statistics_path = (
        "tests/test_data/ascat_output/CNV.somatic.SAMPLE.ascat.samplestatistics.txt"
    )

    # WHEN calling the function
    pdf = add_table_pdf(generate_fpdf(), test_statistics_path)

    # THEN check if the table is appended to the created PDF
    assert isinstance(pdf, FPDF)
    assert pdf.page_no() == 1


def test_create_pdf(tmp_path, cli_runner):
    # GIVEN ascatNGgs output statistics
    statistics_path = (
        "tests/test_data/ascat_output/CNV.somatic.SAMPLE.ascat.samplestatistics.txt"
    )

    # GIVEN ascatNGgs output plots
    plots_path = [
        "tests/test_data/ascat_output/CNV.somatic.SAMPLE.ascat.germline.png",
        "tests/test_data/ascat_output/CNV.somatic.SAMPLE.ascat.sunrise.png",
    ]

    # GIVEN the output path
    output_path = tmp_path / "ascat.output.pdf"

    print(output_path)

    # WHEN invoking the python script
    result = cli_runner.invoke(
        create_pdf, [str(output_path), statistics_path, plots_path[0], plots_path[1]]
    )

    # THEN check if the PDF is correctly created and there is no errors
    assert result.exit_code == 0
    assert Path(output_path).exists()
