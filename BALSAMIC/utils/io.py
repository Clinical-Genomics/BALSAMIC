"""Input/Output utility methods."""
import json
import logging
from datetime import datetime
from pathlib import Path

import snakemake
import yaml
from graphviz import Source

from BALSAMIC import __version__ as balsamic_version
from BALSAMIC.utils.cli import CaptureStdout
from BALSAMIC.utils.exc import BalsamicError

LOG = logging.getLogger(__name__)


def generate_workflow_graph(
    config_path: Path, directory_path: Path, snakefile: Path, title: str
) -> None:
    """Generate snakemake workflow graph and save it in a PDF file."""
    with CaptureStdout() as graph_dot:
        snakemake.snakemake(
            snakefile=snakefile,
            dryrun=True,
            configfiles=[config_path.as_posix()],
            printrulegraph=True,
        )
    graph_title: str = "_".join(["BALSAMIC", balsamic_version, title])
    graph_dot: str = "".join(graph_dot).replace(
        "snakemake_dag {", 'BALSAMIC { label="' + graph_title + '";labelloc="t";'
    )
    graph: Source = Source(
        graph_dot,
        directory=directory_path.as_posix(),
        filename=f"{title}_graph",
        format="pdf",
        engine="dot",
    )
    try:
        graph_pdf: Path = Path(graph.render())
        LOG.info(f"Workflow graph generated successfully ({graph_pdf.as_posix()}) ")
    except Exception:
        LOG.error("Workflow graph generation failed")
        raise BalsamicError()


def read_json(json_path: str) -> dict:
    """Read JSON file and return a dictionary."""
    if Path(json_path).exists():
        with open(json_path, "r") as fn:
            return json.load(fn)
    else:
        raise FileNotFoundError(f"The JSON file {json_path} was not found")


def write_json(json_obj: dict, path: str) -> None:
    """Write JSON format data to an output file."""
    try:
        with open(path, "w") as fn:
            json.dump(json_obj, fn, indent=4)
    except OSError as error:
        raise OSError(f"Error while writing JSON file: {path}, error: {error}")


def read_yaml(yaml_path: str) -> dict:
    """Read data from a yaml file."""
    if Path(yaml_path).exists():
        with open(yaml_path, "r") as fn:
            return yaml.load(fn, Loader=yaml.SafeLoader)
    else:
        raise FileNotFoundError(f"The YAML file {yaml_path} was not found")


def write_finish_file(file_path: str) -> None:
    """Write finish file indicating the analysis completion."""
    with open(file_path, mode="w") as finish_file:
        finish_file.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
