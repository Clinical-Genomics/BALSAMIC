#!/usr/bin/env python
import os
from pathlib import Path
from typing import List, Tuple

import click

from fpdf import FPDF
from PIL import Image


@click.command(short_help="Merge statistics and plots into a single CNV report")
@click.argument("data_paths", nargs=-1, type=click.Path(exists=True), required=False)
@click.option("--output", type=click.Path(exists=False), required=True)
def generate_cnv_report(data_paths: List[Path], output: Path) -> None:
    """Generate a CNV report given a set of statistic files and a list of plots."""
    pdf: PDF = get_pdf_instance()
    statistics, plots = get_pdf_data(data_paths)
    pdf: PDF = add_data_to_pdf(pdf=pdf, data_paths=statistics) if statistics else pdf
    pdf: PDF = add_plots_to_pdf(pdf=pdf, plot_paths=plots) if plots else pdf
    pdf.output(output)


class PDF(FPDF):
    """PDF generation subclass."""

    def footer(self):
        """Overwrite the predetermined method to perform a specific footer processing."""
        self.set_y(-15)
        self.set_font("helvetica", "I", 8)
        self.cell(0, 10, f"Page {self.page_no()}/{{nb}}", 0, 0, "C")


def get_pdf_instance() -> PDF:
    """Return a PDF instance."""
    pdf: PDF = PDF()
    pdf.alias_nb_pages(alias="{nb}")
    return pdf


def get_pdf_data(data_paths: List[Path]) -> Tuple[List[Path], List[Path]]:
    """Return"""
    statistics = []
    plots = []
    for path in data_paths:
        if Path(path).suffix == ".txt":
            statistics.append(path)
        if Path(path).suffix == ".png":
            plots.append(path)
    return statistics, plots


def add_data_to_pdf(pdf: PDF, data_paths: List[Path]) -> PDF:
    """Add statistics to a PDF instance."""
    for data_path in data_paths:
        with open(data_path) as data:
            data = data.readlines()
        pdf.add_page()
        pdf.set_font("helvetica", "B", 15)
        # Title layout & styling
        title: str = os.path.basename(data_path).replace(".txt", "")
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


def add_plots_to_pdf(pdf: PDF, plot_paths: List[Path]) -> PDF:
    """Add plots to a PDF instance."""
    pdf.set_font("helvetica", "B", 15)
    portrait_size = {"sunrise": (500, 500), "circular": (700, 700)}
    portrait_position = {"sunrise": (10, 55), "circular": (-29, 55)}
    for path in plot_paths:
        title: str = os.path.basename(path).replace(".png", "")
        plot_name = title.split(".")[-1]
        # Image & page layout parameters
        if plot_name in portrait_size.keys():
            page_orientation = "portrait"
            img_size = portrait_size[plot_name]
            title_w_pos = 25
            title_wh = 140, 10
            img_xy = portrait_position[plot_name]
        else:
            page_orientation = "landscape"
            img_size = (800, 800) if "scatter" != plot_name else (650, 650)
            title_w_pos = 68.5
            title_wh = 140, 10
            img_xy = (5, 40) if "scatter" != plot_name else (25, 23)

        pdf.add_page(orientation=page_orientation)

        # Title position & styling
        pdf.cell(title_w_pos)
        pdf.cell(title_wh[0], title_wh[1], title, 1, 0, "C")

        # Image position & resizing
        img = Image.open(path)
        img.thumbnail(img_size, Image.ANTIALIAS)
        pdf.image(img, img_xy[0], img_xy[1])

    return pdf


if __name__ == "__main__":
    generate_cnv_report()
