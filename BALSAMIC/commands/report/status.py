import os
import logging
import json
import click
import snakemake
from colorclass import Color

from BALSAMIC.utils.cli import get_snakefile
from BALSAMIC.utils.cli import CaptureStdout
from BALSAMIC.utils.cli import get_file_status_string
from BALSAMIC.utils.rule import get_result_dir
from BALSAMIC.utils.constants import VCF_DICT

LOG = logging.getLogger(__name__)


@click.command(
    "status",
    short_help=
    "Creates a YAML file with output from variant caller and alignment.",
)
@click.option(
    "-s",
    "--sample-config",
    required=True,
    help="Sample config file. Output of balsamic config sample",
)
@click.option("-m",
              "--show-only-missing",
              is_flag=True,
              default=False,
              show_default=True,
              help="Only show missing files.")
@click.option(
    "-p",
    "--print-files",
    is_flag=True,
    default=False,
    show_default=True,
    help="Print list of files. Otherwise only final count will be printed.")
@click.option(
    '--disable-variant-caller',
    help=
    f'Run workflow with selected variant caller(s) disable. Use comma to remove multiple variant callers. Valid '
    f'values are: {list(VCF_DICT.keys())}',
)
@click.pass_context
def status(context, sample_config, show_only_missing, print_files,
           disable_variant_caller):
    """
    cli for status sub-command.
    """
    LOG.info(f"BALSAMIC started with log level {context.obj['loglevel']}.")
    LOG.debug("Reading input sample config")
    with open(sample_config, "r") as fn:
        sample_config_dict = json.load(fn)

    result_dir = get_result_dir(sample_config_dict)
    analysis_type = sample_config_dict["analysis"]["analysis_type"]
    sequencing_type = sample_config_dict["analysis"]["sequencing_type"]
    snakefile = get_snakefile(analysis_type, sequencing_type)

    with CaptureStdout() as summary:
        snakemake.snakemake(
            snakefile=snakefile,
            dryrun=True,
            summary=True,
            config={"disable_variant_caller": disable_variant_caller},
            configfiles=[sample_config],
            quiet=True,
        )
    summary = [i.split("\t") for i in summary]
    summary_dict = [dict(zip(summary[0], value)) for value in summary[1:]]

    if not os.path.isfile(os.path.join(result_dir, "analysis_finish")):
        LOG.warning(
            "analysis_finish file is missing. Analysis might be incomplete or running."
        )

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

    finish_file_count = 'Finished file count: {}'.format(len(existing_files))
    missing_file_count = 'Missing file count: {}'.format(len(missing_files))
    click.echo(Color('{yellow}Final tally:{/yellow}'))
    click.echo(Color('{yellow}\t' + finish_file_count + '{/yellow}'))
    click.echo(Color('{yellow}\t' + missing_file_count + '{/yellow}'))
