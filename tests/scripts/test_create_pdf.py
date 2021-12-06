from fpdf import FPDF

from BALSAMIC.assets.scripts.create_pdf import (
    generate_fpdf,
    add_images_pdf,
    add_table_pdf,
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
