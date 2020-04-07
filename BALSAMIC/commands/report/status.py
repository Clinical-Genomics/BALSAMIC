import os
import sys
import logging
import glob
import json
import yaml
import click
import copy
import snakemake
from collections import defaultdict
from yapf.yapflib.yapf_api import FormatFile

from BALSAMIC.utils.cli import get_from_two_key
from BALSAMIC.utils.cli import merge_dict_on_key
from BALSAMIC.utils.cli import get_file_extension
from BALSAMIC.utils.cli import find_file_index
from BALSAMIC.utils.cli import write_json
from BALSAMIC.utils.cli import get_snakefile
from BALSAMIC.utils.cli import CaptureStdout
from BALSAMIC.utils.cli import get_file_status_string
from BALSAMIC.utils.rule import get_result_dir
from BALSAMIC.utils.exc import BalsamicError

LOG = logging.getLogger(__name__)


@click.command(
    "status",
    short_help="Creates a YAML file with output from variant caller and alignment.",
)
@click.option(
    "--sample-config",
    required=True,
    help="Sample config file. Output of balsamic config sample",
)
@click.option("--only-missing", is_flag=True,
    show_default=True,
    help="Only show missing files.")
@click.pass_context
def status(context, sample_config, only_missing):
    """
    cli for deliver sub-command.
    Writes <case_id>.hk in result_directory.
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
            config={"delivery": "True"},
            dryrun=True,
            summary=True,
            configfiles=[sample_config],
            quiet=True,
        )
    summary = [i.split("\t") for i in summary]
    summary_dict = [dict(zip(summary[0], value)) for value in summary[1:]]

    if not os.path.isfile(os.path.join(result_dir, "analysis_finish")):
        LOG.warning("analysis_finish file is missing. Analysis might be incomplete") 

        
    for entries in summary_dict:
        delivery_file = entries["output_file"]
  
        file_stat_str = get_file_status_string(delivery_file, only_missing) 
        if file_stat_str:
            click.echo(file_stat_str)
