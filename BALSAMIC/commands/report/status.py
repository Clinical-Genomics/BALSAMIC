"""Balsamic status report CLI."""
import json
import logging
import os

import click
import snakemake
from colorclass import Color

from BALSAMIC.commands.options import OPTION_SAMPLE_CONFIG
from BALSAMIC.commands.report.options import (
    OPTION_PRINT_FILES,
    OPTION_SHOW_ONLY_MISSING_FILES,
)
from BALSAMIC.utils.cli import CaptureStdout
from BALSAMIC.utils.cli import get_file_status_string
from BALSAMIC.utils.cli import get_snakefile
from BALSAMIC.utils.rule import get_result_dir

LOG = logging.getLogger(__name__)


@click.command("status", short_help="Print the analysis file status.")
@OPTION_SAMPLE_CONFIG
@OPTION_SHOW_ONLY_MISSING_FILES
@OPTION_PRINT_FILES
@click.pass_context
def status(
    context: click.Context,
    sample_config: str,
    show_only_missing: bool,
    print_files: bool,
):
    """Analysis status CLI command."""
    LOG.info(f"BALSAMIC started with log level {context.obj['log_level']}.")
    LOG.debug("Reading input sample config")
    with open(sample_config, "r") as fn:
        sample_config_dict = json.load(fn)

    result_dir = get_result_dir(sample_config_dict)
    analysis_type = sample_config_dict["analysis"]["analysis_type"]
    analysis_workflow = sample_config_dict["analysis"]["analysis_workflow"]
    snakefile = get_snakefile(analysis_type, analysis_workflow)

    if os.path.isfile(os.path.join(result_dir, "analysis_finish")):
        snakemake.snakemake(
            snakefile=snakefile,
            config={
                "benchmark_plots": "True",
            },
            dryrun=True,
            configfiles=[sample_config],
            quiet=True,
        )
    else:
        LOG.warning(
            "analysis_finish file is missing. Analysis might be incomplete or running."
        )

    with CaptureStdout() as summary:
        snakemake.snakemake(
            snakefile=snakefile,
            dryrun=True,
            summary=True,
            configfiles=[sample_config],
            quiet=True,
        )
    summary = [i.split("\t") for i in summary]
    summary_dict = [dict(zip(summary[1], value)) for value in summary[2:]]

    existing_files = set()
    missing_files = set()

    for entries in summary_dict:
        delivery_file = entries["output_file"]

        file_status_str, file_status = get_file_status_string(delivery_file)
        if file_status and print_files:
            click.echo(file_status_str)

        if not file_status and (show_only_missing or print_files):
            click.echo(file_status_str)

        if file_status:
            existing_files.add(delivery_file)
        if not file_status:
            missing_files.add(delivery_file)

    finish_file_count = "Finished file count: {}".format(len(existing_files))
    missing_file_count = "Missing file count: {}".format(len(missing_files))
    click.echo(Color("{yellow}Final tally:{/yellow}"))
    click.echo(Color("{yellow}\t" + finish_file_count + "{/yellow}"))
    click.echo(Color("{yellow}\t" + missing_file_count + "{/yellow}"))
