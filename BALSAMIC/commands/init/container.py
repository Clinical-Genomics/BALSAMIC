import os
import logging
import subprocess
from pathlib import Path

import click
import graphviz
import snakemake

from BALSAMIC.utils.cli import write_json, merge_json, CaptureStdout, get_snakefile, SnakeMake
from BALSAMIC import __version__ as balsamic_version

LOG = logging.getLogger(__name__)


@click.command("container",
               short_help="Download matching version for container")
@click.option("-o",
              "--outdir",
              "--out-dir",
              required=True,
              help="Output directory for container files.")
@click.option("-v",
              "--balsamic-version",
              default=balsamic_version,
              help="Container for BALSAMIC version to download")
@click.option('-f',
              '--force',
              show_default=True,
              default=False,
              is_flag=True,
              help="Force re-downloading all containers")
@click.pass_context
def container(context, outdir, force):
    pass
