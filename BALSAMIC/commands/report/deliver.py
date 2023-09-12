"""Balsamic report delivery CLI."""
import json
import logging
import os
import subprocess
import sys
from pathlib import Path
from typing import List

import click
import snakemake
import yaml

from BALSAMIC.commands.options import (
    OPTION_SAMPLE_CONFIG,
    OPTION_DISABLE_VARIANT_CALLER,
    OPTION_DELIVERY_MODE,
    OPTION_RULES_TO_DELIVER,
)
from BALSAMIC.constants.analysis import RunMode, RuleDeliveryMode
from BALSAMIC.constants.rules import DELIVERY_RULES
from BALSAMIC.models.snakemake import SnakemakeExecutable
from BALSAMIC.utils.cli import convert_deliverables_tags
from BALSAMIC.utils.cli import get_file_extension
from BALSAMIC.utils.cli import get_snakefile
from BALSAMIC.utils.io import write_json
from BALSAMIC.utils.rule import get_result_dir

LOG = logging.getLogger(__name__)


@click.command("deliver", short_help="Creates a report file with output files")
@OPTION_DELIVERY_MODE
@OPTION_DISABLE_VARIANT_CALLER
@OPTION_RULES_TO_DELIVER
@OPTION_SAMPLE_CONFIG
@click.pass_context
def deliver(
    context: click.Context,
    delivery_mode: RuleDeliveryMode,
    disable_variant_caller: str,
    rules_to_deliver: List[str],
    sample_config: str,
):
    """Deliver command to write <case_id>.hk with the output analysis files."""
    LOG.info(f"BALSAMIC started with log level {context.obj['log_level']}.")
    LOG.debug("Reading input sample config")
    with open(sample_config, "r") as fn:
        sample_config_dict = json.load(fn)

    default_rules_to_deliver = DELIVERY_RULES

    if not rules_to_deliver:
        rules_to_deliver = default_rules_to_deliver

    rules_to_deliver = list(rules_to_deliver)
    if delivery_mode == RuleDeliveryMode.APPEND:
        rules_to_deliver.extend(default_rules_to_deliver)

    case_name = sample_config_dict["analysis"]["case_id"]
    result_dir = get_result_dir(sample_config_dict)
    dst_directory = os.path.join(result_dir, "delivery_report")
    LOG.info("Creating delivery_report directory")
    Path.mkdir(Path(dst_directory), parents=True, exist_ok=True)

    yaml_write_directory = os.path.join(result_dir, "delivery_report")
    Path.mkdir(Path(yaml_write_directory), parents=True, exist_ok=True)

    analysis_type = sample_config_dict["analysis"]["analysis_type"]
    analysis_workflow = sample_config_dict["analysis"]["analysis_workflow"]
    snakefile = get_snakefile(analysis_type, analysis_workflow)

    report_path = Path(yaml_write_directory, f"{case_name}_report.html")
    LOG.info(f"Creating report file {report_path.as_posix()}")

    LOG.info(f"Delivering {analysis_workflow} workflow...")
    working_dir = Path(
        sample_config_dict["analysis"]["analysis_dir"], case_name, "BALSAMIC_run"
    )
    snakemake_executable: SnakemakeExecutable = SnakemakeExecutable(
        case_id=case_name,
        config_path=sample_config,
        disable_variant_caller=disable_variant_caller,
        report_path=report_path,
        run_analysis=True,
        run_mode=RunMode.LOCAL,
        snakefile=snakefile,
        snakemake_options=["--quiet"],
        working_dir=working_dir,
    )
    subprocess.check_output(
        f"{sys.executable} -m {snakemake_executable.get_command()}".split(),
        shell=False,
    )
    LOG.info(f"Workflow report file {report_path.as_posix()}")

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
            "path": report_path.as_posix(),
            "step": "balsamic_delivery",
            "format": get_file_extension(report_path.as_posix()),
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
