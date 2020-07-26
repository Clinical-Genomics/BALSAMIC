import os
import sys
import logging
import glob
import json
import yaml
import click
import copy
import snakemake
import datetime
import subprocess
from collections import defaultdict
from yapf.yapflib.yapf_api import FormatFile
from pathlib import Path

from BALSAMIC.utils.cli import get_from_two_key
from BALSAMIC.utils.cli import merge_dict_on_key
from BALSAMIC.utils.cli import get_file_extension
from BALSAMIC.utils.cli import find_file_index
from BALSAMIC.utils.cli import write_json
from BALSAMIC.utils.cli import get_snakefile
from BALSAMIC.utils.cli import CaptureStdout
from BALSAMIC.utils.cli import SnakeMake
from BALSAMIC.utils.rule import get_result_dir
from BALSAMIC.utils.exc import BalsamicError

LOG = logging.getLogger(__name__)


@click.command(
    "deliver",
    short_help=
    "Creates a YAML file with output from variant caller and alignment.",
)
@click.option(
    "--sample-config",
    "-s",
    required=True,
    help="Sample config file. Output of balsamic config sample",
)
@click.option('-a',
              '--analysis-type',
              required=False,
              type=click.Choice(['qc', 'paired', 'single']),
              help='Type of analysis to run from input config file.\
              By default it will read from config file, but it will override config file \
              if it is set here.')
@click.pass_context
def deliver(context, sample_config, analysis_type):
    """
    cli for deliver sub-command.
    Writes <case_id>.hk in result_directory.
    """
    LOG.info(f"BALSAMIC started with log level {context.obj['loglevel']}.")
    LOG.debug("Reading input sample config")
    with open(sample_config, "r") as fn:
        sample_config_dict = json.load(fn)

    case_name = sample_config_dict['analysis']['case_id']
    result_dir = get_result_dir(sample_config_dict)
    dst_directory = os.path.join(result_dir, "delivery_report")
    LOG.info("Creatiing delivery_report directory")
    os.makedirs(dst_directory, exist_ok=True)

    yaml_write_directory = os.path.join(result_dir, "delivery_report")
    os.makedirs(yaml_write_directory, exist_ok=True)

    analysis_type = analysis_type if analysis_type else sample_config_dict['analysis']['analysis_type']
    sequencing_type = sample_config_dict["analysis"]["sequencing_type"]
    snakefile = get_snakefile(analysis_type, sequencing_type)

    report_file_name = os.path.join(
        yaml_write_directory, sample_config_dict["analysis"]["case_id"] + "_report.html"
    )
    LOG.info("Creating report file {}".format(report_file_name))

    # write report.html file
    report = SnakeMake()
    report.case_name = case_name 
    report.working_dir = os.path.join(sample_config_dict['analysis']['analysis_dir'] , \
        sample_config_dict['analysis']['case_id'], 'BALSAMIC_run')
    report.report = report_file_name
    report.configfile = sample_config
    report.snakefile = snakefile 
    report.run_mode = 'local'
    report.use_singularity = False
    report.run_analysis = True
    cmd=sys.executable + " -m  " + report.build_cmd()
    subprocess.check_output(cmd.split(), shell=False)

    snakemake.snakemake(
        snakefile=snakefile,
        config={"delivery": "True", "rules_to_deliver": "multiqc,vep_somatic,vep_germline,vep_stat,ngs_filter_vardict"},
        dryrun=True,
        configfiles=[sample_config],
        quiet=True,
    )

    delivery_file_name = os.path.join(
        yaml_write_directory,
        case_name + ".hk")

    delivery_file_ready = os.path.join(
        yaml_write_directory,
        case_name + "_delivery_ready.hk",
    )
    with open(delivery_file_ready, "r") as fn:
        delivery_file_ready_dict = json.load(fn)

    delivery_json = dict()
    delivery_json["files"] = list()

    # Add Housekeeper file to report
    delivery_json["files"].append({
        "path":
        report_file_name,
        "step":
        "balsamic_delivery",
        "format":
        "html",
        "tag":
        "report",
        "id":
        case_name,
    })
    # Add CASE_ID.JSON to report
    delivery_json["files"].append({
        "path":
        Path(sample_config).resolve().as_posix(),
        "step":
        "case_config",
        "format":
        "json",
        "tag":
        "config",
        "id":
        case_name,
    })
    # Add DAG Graph to report
    delivery_json["files"].append({
        "path":
        sample_config_dict["analysis"]["dag"],
        "step":
        "case_config",
        "format":
        "pdf",
        "tag":
        "dag",
        "id":
        case_name,
    })

    delivery_json["files"].extend(delivery_file_ready_dict)

    write_json(delivery_json, delivery_file_name)
    with open(delivery_file_name + ".yaml", "w") as fn:
        yaml.dump(delivery_json, fn, default_flow_style=False)

    LOG.info(f"Housekeeper delivery file {delivery_file_name}")
    LOG.info(f"Workflow report file {report_file_name}")
