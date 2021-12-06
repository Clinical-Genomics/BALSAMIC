import os
import subprocess
import json
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import h5py
import typing

from fpdf import FPDF
from PIL import Image

from BALSAMIC.utils.rule import get_threads
from BALSAMIC.utils.cli import get_config
from BALSAMIC.utils.cli import generate_h5


def plot_analysis(
        log_file: Path, h5_file: Path, fig_name: Path
) -> typing.Union[None, Path]:
    """
    plots analysis job.
    """

    cluster_config = get_config("cluster")
    with open(cluster_config, "r") as f:
        cluster_config = json.load(f)

    log_file_list = Path(log_file).name.split(".")

    job_name = ".".join(log_file_list[0:4])
    rule_name = log_file_list[2]
    mem_per_core = 5222
    requested_cores = get_threads(cluster_config, rule_name)
    case_name = log_file_list[1]
    job_id = log_file_list[4].split("_")[1]

    # This is lazy and memory inefficient, but it gets the job done.
    df_array = h5py.File(h5_file, "r")
    node_name = list(df_array["Steps"]["batch"]["Nodes"].keys())[0]

    if not "Tasks" in list(df_array["Steps"]["batch"]["Nodes"][node_name]):
        return None

    df = pd.DataFrame(
        np.array(df_array["Steps"]["batch"]["Nodes"][node_name]["Tasks"]["0"])
    )

    # Convert kilohurtz to gigahurtz
    df["CPUFrequency"] = df["CPUFrequency"] / 1e6

    # Convert kb to gb
    df["RSS"] = df["RSS"] / 1e6
    df["VMSize"] = df["VMSize"] / 1e6

    figure_title = "Case name: {}\nRule: {}\nRun time: {} seconds\nJob name: {}\nJob ID: {}".format(
        case_name, rule_name, df["ElapsedTime"].iloc[-1], job_name, job_id
    )

    plt.rcParams["figure.figsize"] = [10, 10]

    fig, (cpu_ax, mem_ax, io_ax) = plt.subplots(nrows=3)
    fig.suptitle(figure_title, fontsize=12, horizontalalignment="center")

    cpu_ax_color = "b"
    df.plot(
        y="CPUUtilization", x="ElapsedTime", ax=cpu_ax, color=cpu_ax_color, style="--"
    )
    cpu_ax.set_title("CPU statistics")
    cpu_ax.set_xlabel("Wall seconds")
    cpu_ax.set_ylabel("Core usage (max {}%)".format(requested_cores * 100))
    cpu_ax.yaxis.label.set_color(cpu_ax_color)
    cpu_ax.yaxis.label.set_color(cpu_ax_color)
    cpu_ax.tick_params(axis="y", colors=cpu_ax_color)
    cpu_ax.legend(loc="best", frameon=False)
    cpu_ax.spines["top"].set_visible(False)
    cpu_ax.spines["right"].set_visible(False)
    max_cpu_line = cpu_ax.axhline(requested_cores * 100, color=cpu_ax_color, ls="-")
    max_cpu_line.set_label("Max available")

    mem_ax_color = "g"
    df.plot(y="VMSize", x="ElapsedTime", ax=mem_ax, color=mem_ax_color, style="--")
    mem_ax.set_title("Memory statistics")
    mem_ax.set_xlabel("Wall seconds")
    mem_ax.set_ylabel(
        "Memory usage GB (max {}GB)".format(
            round(mem_per_core * requested_cores / 1024)
        )
    )
    mem_ax.yaxis.label.set_color(mem_ax_color)
    mem_ax.yaxis.label.set_color(mem_ax_color)
    mem_ax.tick_params(axis="y", colors=mem_ax_color)
    mem_ax.legend(loc="best", frameon=False)
    mem_ax.spines["top"].set_visible(False)
    mem_ax.spines["right"].set_visible(False)
    max_cpu_line = mem_ax.axhline(
        round(mem_per_core * requested_cores / 1024), color=mem_ax_color, ls="-"
    )
    max_cpu_line.set_label("Max available mem")

    read_io_ax_color = "m"
    read_io_ax = df.plot(
        y="ReadMB",
        x="ElapsedTime",
        color=read_io_ax_color,
        style="--",
        ax=io_ax,
        legend=False,
    )
    read_io_ax.set_xlabel("Wall seconds")
    read_io_ax.set_ylabel("Disk read (MANIFEST.inB)")
    read_io_ax.yaxis.label.set_color(read_io_ax_color)
    read_io_ax.yaxis.label.set_color(read_io_ax_color)
    read_io_ax.tick_params(axis="y", colors=read_io_ax_color)
    read_io_ax.spines["top"].set_visible(False)

    write_io_ax = read_io_ax.twinx()
    write_io_ax_color = "olive"
    write_io_ax = df.plot(
        y="WriteMB",
        x="ElapsedTime",
        ax=write_io_ax,
        color=write_io_ax_color,
        style="--",
        legend=False,
    )
    write_io_ax.set_title("Disk I/O statistics")
    write_io_ax.set_xlabel("Wall seconds")
    write_io_ax.set_ylabel("Disk write (MB)")
    write_io_ax.yaxis.label.set_color(write_io_ax_color)
    write_io_ax.yaxis.label.set_color(write_io_ax_color)
    write_io_ax.tick_params(axis="y", colors=write_io_ax_color)
    write_io_ax.yaxis.tick_right()
    write_io_ax.spines["top"].set_visible(False)

    handles, labels = [], []
    for ax in [write_io_ax, read_io_ax]:
        for h, l in zip(*ax.get_legend_handles_labels()):
            handles.append(h)
            labels.append(l)

    plt.legend(handles, labels, loc="best", ncol=len(handles), frameon=False)

    plt.tight_layout()
    plt.savefig(fig_name, dpi=300)
    plt.close()
    return fig_name


class PDF(FPDF):
    def footer(self):
        self.set_y(-15)
        self.set_font("helvetica", "I", 8)
        self.cell(0, 10, f"Page {self.page_no()}/{{nb}}", 0, 0, "C")


def create_pdf():
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


def save_ascat_output_pdf(output_path, data_path, img_paths):
    pdf = create_pdf()
    pdf = add_table_pdf(pdf, data_path)
    pdf = add_images_pdf(pdf, img_paths)
    pdf.output(output_path)
