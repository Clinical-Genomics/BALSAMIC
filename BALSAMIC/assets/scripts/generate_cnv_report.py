#!/usr/bin/env python
import os
import click

from fpdf import FPDF
from PIL import Image


@click.command(
    short_help="Merge images and a txt file into a single PDF",
)
@click.argument("output", type=click.Path(exists=False), required=True)
@click.argument("data", type=click.Path(exists=True), required=True)
@click.argument("images", nargs=-1, type=click.Path(exists=True), required=True)
def generate_cnv_report(output, data, images):
    pdf = generate_fpdf()
    pdf = add_table_pdf(pdf, data)
    pdf = add_images_pdf(pdf, images)
    pdf.output(output)


class PDF(FPDF):
    def footer(self):
        self.set_y(-15)
        self.set_font("helvetica", "I", 8)
        self.cell(0, 10, f"Page {self.page_no()}/{{nb}}", 0, 0, "C")


def generate_fpdf():
    pdf = PDF()
    pdf.alias_nb_pages(alias="{nb}")
    return pdf


def add_images_pdf(pdf, img_paths):
    pdf.set_font("helvetica", "B", 15)

    for path in img_paths:
        title = os.path.basename(path).replace(".png", "")

        # Image & page layout parameters
        if "sunrise" in title:
            page_orientation = "portrait"
            img_size = 500, 500
            title_w_pos = 25
            title_wh = 140, 10
            img_xy = 10, 55
        else:
            page_orientation = "landscape"
            img_size = 800, 800
            title_w_pos = 68.5
            title_wh = 140, 10
            img_xy = 5, 40

        pdf.add_page(orientation=page_orientation)

        # Title position & styling
        pdf.cell(title_w_pos)
        pdf.cell(title_wh[0], title_wh[1], title, 1, 0, "C")

        # Image position & resizing
        img = Image.open(path)
        img.thumbnail(img_size, Image.ANTIALIAS)
        pdf.image(img, img_xy[0], img_xy[1])

    return pdf


def add_table_pdf(pdf, data_path):

    with open(data_path) as data:
        data = data.readlines()

    pdf.add_page()
    pdf.set_font("helvetica", "B", 15)

    # Title layout & styling
    title = os.path.basename(data_path).replace(".txt", "")
    pdf.cell(25)
    pdf.cell(140, 10, title, 1, 0, "C")
    pdf.cell(35, 25, ln=1)  # Post title indentation

    # Table layout & styling
    pdf.set_font("Times", size=11)
    line_height = pdf.font_size * 2.5
    col_width = pdf.epw / 4  # Even distribution of the content
    for row in data:
        pdf.cell(45)
        for statistic in row.split():
            pdf.multi_cell(
                col_width,
                line_height,
                statistic,
                align="C",
                border=1,
                ln=3,
                max_line_height=pdf.font_size,
            )
        pdf.ln(line_height)

    return pdf


if __name__ == "__main__":
    generate_cnv_report()
