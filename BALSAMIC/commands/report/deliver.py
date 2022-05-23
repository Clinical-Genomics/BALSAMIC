import os
import sys
import logging
import json
import yaml
import click
import snakemake
import subprocess
from pathlib import Path

from BALSAMIC.utils.cli import get_file_extension, read_yaml
from BALSAMIC.utils.cli import write_json
from BALSAMIC.utils.cli import get_snakefile
from BALSAMIC.utils.cli import SnakeMake
from BALSAMIC.utils.cli import convert_deliverables_tags
from BALSAMIC.utils.rule import get_result_dir
from BALSAMIC.utils.exc import BalsamicError
from BALSAMIC.constants.workflow_params import VCF_DICT
from BALSAMIC.constants.workflow_rules import DELIVERY_RULES

LOG = logging.getLogger(__name__)


@click.command(
    "deliver",
    short_help="Creates a YAML file with output from variant caller and alignment.",
)
@click.option(
    "--sample-config",
    "-s",
    required=True,
    help="Sample config file. Output of balsamic config sample",
)
@click.option(
    "-a",
    "--analysis-type",
    required=False,
    type=click.Choice(["paired", "single"]),
    help=(
        "Type of analysis to run from input config file."
        "By default it will read from config file, but it will override config file"
        "if it is set here."
    ),
)
@click.option(
    "-r",
    "--rules-to-deliver",
    multiple=True,
    help=f"Specify a rule to deliver. Delivery mode selected via --delivery-mode option."
    f"Current available rules to deliver are: {', '.join(DELIVERY_RULES)} ",
)
@click.option(
    "-m",
    "--delivery-mode",
    type=click.Choice(["a", "r"]),
    default="a",
    show_default=True,
    help="a: append rules-to-deliver to current delivery options. "
    "r: reset current rules to delivery to only the ones specified",
)
@click.option(
    "--disable-variant-caller",
    help=f"Run workflow with selected variant caller(s) disable. Use comma to remove multiple variant callers. Valid "
    f"values are: {list(VCF_DICT.keys())}",
)
@click.pass_context
def deliver(
    context,
    sample_config,
    analysis_type,
    rules_to_deliver,
    delivery_mode,
    disable_variant_caller,
):
    """
    cli for deliver sub-command.
    Writes <case_id>.hk in result_directory.
    """
    LOG.info(f"BALSAMIC started with log level {context.obj['loglevel']}.")
    LOG.debug("Reading input sample config")
    with open(sample_config, "r") as fn:
        sample_config_dict = json.load(fn)

    default_rules_to_deliver = DELIVERY_RULES

    if not rules_to_deliver:
        rules_to_deliver = default_rules_to_deliver

    rules_to_deliver = list(rules_to_deliver)
    if delivery_mode == "a":
        rules_to_deliver.extend(default_rules_to_deliver)

    case_name = sample_config_dict["analysis"]["case_id"]
    result_dir = get_result_dir(sample_config_dict)
    dst_directory = os.path.join(result_dir, "delivery_report")
    LOG.info("Creating delivery_report directory")
    Path.mkdir(Path(dst_directory), parents=True, exist_ok=True)

    yaml_write_directory = os.path.join(result_dir, "delivery_report")
    Path.mkdir(Path(yaml_write_directory), parents=True, exist_ok=True)

    analysis_type = (
        analysis_type
        if analysis_type
        else sample_config_dict["analysis"]["analysis_type"]
    )
    reference_genome = sample_config_dict["reference"]["reference_genome"]
    analysis_workflow = sample_config_dict["analysis"]["analysis_workflow"]
    snakefile = get_snakefile(analysis_type, analysis_workflow, reference_genome)

    report_file_name = os.path.join(
        yaml_write_directory, sample_config_dict["analysis"]["case_id"] + "_report.html"
    )
    LOG.info("Creating report file {}".format(report_file_name))

    # write report.html file
    report = SnakeMake()
    report.case_name = case_name
    report.working_dir = os.path.join(
        sample_config_dict["analysis"]["analysis_dir"],
        sample_config_dict["analysis"]["case_id"],
        "BALSAMIC_run",
    )
    report.report = report_file_name
    report.configfile = sample_config
    report.snakefile = snakefile
    report.run_mode = "local"
    report.use_singularity = False
    report.run_analysis = True
    report.sm_opt = ["--quiet"]
    if disable_variant_caller:
        report.disable_variant_caller = disable_variant_caller
    cmd = sys.executable + " -m  " + report.build_cmd()
    subprocess.check_output(cmd.split(), shell=False)
    LOG.info(f"Workflow report file {report_file_name}")

    snakemake.snakemake(
        snakefile=snakefile,
        config={"delivery": "True", "rules_to_deliver": ",".join(rules_to_deliver)},
        dryrun=True,
        configfiles=[sample_config],
        quiet=True,
    )

    delivery_file_name = os.path.join(yaml_write_directory, case_name + ".hk")

    delivery_file_ready = os.path.join(
        yaml_write_directory,
        case_name + "_delivery_ready.hk",
    )
    with open(delivery_file_ready, "r") as fn:
        delivery_file_ready_dict = json.load(fn)

    delivery_json = dict()
    delivery_json["files"] = delivery_file_ready_dict

    delivery_json = convert_deliverables_tags(
        delivery_json=delivery_json, sample_config_dict=sample_config_dict
    )

    # Add Housekeeper file to report
    delivery_json["files"].append(
        {
            "path": report_file_name,
            "step": "balsamic_delivery",
            "format": get_file_extension(report_file_name),
            "tag": ["balsamic-report"],
            "id": case_name,
        }
    )
    # Add CASE_ID.JSON to report
    delivery_json["files"].append(
        {
            "path": Path(sample_config).resolve().as_posix(),
            "step": "case_config",
            "format": get_file_extension(sample_config),
            "tag": ["balsamic-config"],
            "id": case_name,
        }
    )
    # Add DAG Graph to report
    delivery_json["files"].append(
        {
            "path": sample_config_dict["analysis"]["dag"],
            "step": "case_config",
            "format": get_file_extension(sample_config_dict["analysis"]["dag"]),
            "tag": ["balsamic-dag"],
            "id": case_name,
        }
    )

    write_json(delivery_json, delivery_file_name)
    with open(delivery_file_name + ".yaml", "w") as fn:
        yaml.dump(delivery_json, fn, default_flow_style=False)

    LOG.info(f"Housekeeper delivery file {delivery_file_name}")
