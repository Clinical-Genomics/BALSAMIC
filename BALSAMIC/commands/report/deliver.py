import os
import sys
import logging
import json
import yaml
import click
import snakemake
import subprocess
from pathlib import Path

from BALSAMIC.utils.cli import get_file_extension
from BALSAMIC.utils.cli import write_json
from BALSAMIC.utils.cli import get_snakefile
from BALSAMIC.utils.cli import SnakeMake
from BALSAMIC.utils.cli import convert_deliverables_tags
from BALSAMIC.utils.rule import get_result_dir
from BALSAMIC.utils.constants import VCF_DICT
from BALSAMIC.utils.exc import BalsamicError
from BALSAMIC.utils.qc_metrics import get_qc_metrics
from BALSAMIC.utils.qc_report import render_html, report_data_population

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
@click.option(
    "--sample-id-map",
    required=True,
    help=("Separated internal sample ID with external ID. Use comma for"
          "multiple samples. These IDs MUST exist in sample-config."
          "Syntax: internal_id:sample_type:external_id"
          ". e.g. ACC1:tumor:KS454,ACC2:normal:KS556"))
@click.option("--case-id-map",
              required=True,
              help=("Separated internal case ID with external ID."
                    "Syntax: gene_panel_name:external_id"
                    ". e.g. gmck-solid:KSK899:apptag"))
@click.option(
    "-a",
    "--analysis-type",
    required=False,
    type=click.Choice(["qc", "paired", "single"]),
    help=(
        "Type of analysis to run from input config file."
        "By default it will read from config file, but it will override config file"
        "if it is set here."),
)
@click.option(
    "-r",
    "--rules-to-deliver",
    multiple=True,
    help=("Specify a rule to deliver. Delivery "
          "mode selected via --delivery-mode option"),
)
@click.option(
    "-m",
    "--delivery-mode",
    type=click.Choice(["a", "r"]),
    default="a",
    show_default=True,
    help=(
        'a: append rules-to-deliver to current delivery '
        'options. or r: reset current rules to delivery to only the ones specified'
    ))
@click.option(
    '--disable-variant-caller',
    help=
    f'Run workflow with selected variant caller(s) disable. Use comma to remove multiple variant callers. Valid '
    f'values are: {list(VCF_DICT.keys())}',
)
@click.pass_context
def deliver(context, sample_config, analysis_type, rules_to_deliver,
            delivery_mode, disable_variant_caller, sample_id_map, case_id_map):
    """
    cli for deliver sub-command.
    Writes <case_id>.hk in result_directory.
    """
    LOG.info(f"BALSAMIC started with log level {context.obj['loglevel']}.")
    LOG.debug("Reading input sample config")
    with open(sample_config, "r") as fn:
        sample_config_dict = json.load(fn)

    default_rules_to_deliver = [
        "fastp",
        "multiqc",
        "vep_somatic",
        "vep_germline",
        "vep_stat",
        "ngs_filter_vardict",
        "mergeBam_tumor",
        "mergeBam_normal",
        "cnvkit_paired",
        "cnvkit_single",
        "sentieon_dedup",
    ]

    if not rules_to_deliver:
        rules_to_deliver = default_rules_to_deliver

    rules_to_deliver = list(rules_to_deliver)
    if delivery_mode == "a":
        rules_to_deliver.extend(default_rules_to_deliver)

    case_name = sample_config_dict["analysis"]["case_id"]
    result_dir = get_result_dir(sample_config_dict)
    dst_directory = os.path.join(result_dir, "delivery_report")
    LOG.info("Creatiing delivery_report directory")
    os.makedirs(dst_directory, exist_ok=True)

    yaml_write_directory = os.path.join(result_dir, "delivery_report")
    os.makedirs(yaml_write_directory, exist_ok=True)

    analysis_type = (analysis_type if analysis_type else
                     sample_config_dict["analysis"]["analysis_type"])
    sequencing_type = sample_config_dict["analysis"]["sequencing_type"]
    snakefile = get_snakefile(analysis_type, sequencing_type)

    if sequencing_type != "wgs":
        case_id_map = case_id_map.split(":")
        sample_id_map = sample_id_map.split(",")
        sample_map = dict()
        sample_type = dict()
        for sample in sample_id_map:
            lims_id = sample.split(":")[0]
            sample_map[lims_id] = sample.split(":")[1]
            sample_type[lims_id] = sample.split(":")[2]

        meta = dict()
        meta["sample_map"] = sample_map
        meta["sample_type"] = sample_type
        meta["now"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
        meta["config_date"] = sample_config_dict["analysis"][
            "config_creation_date"]
        meta["internal_case_id"] = case_name
        meta["gene_panel_name"] = case_id_map[0]
        meta["case_name"] = case_id_map[1]
        meta["apptag"] = case_id_map[2]

        collected_qc = get_qc_metrics(sample_config_dict["analysis"]["result"])
        meta = report_data_population(collected_qc=collected_qc, meta=meta)
        balsamic_qc_report = os.path.join(yaml_write_directory,
                                          case_name + "_qc_report.html")
        balsamic_qc_report = render_html(meta=meta,
                                         html_out=balsamic_qc_report)

    report_file_name = os.path.join(
        yaml_write_directory,
        sample_config_dict["analysis"]["case_id"] + "_report.html")
    LOG.info("Creating report file {}".format(report_file_name))

    # write report.html file
    report = SnakeMake()
    report.case_name = case_name
    report.working_dir = os.path.join(
        sample_config_dict['analysis']['analysis_dir'],
        sample_config_dict['analysis']['case_id'], 'BALSAMIC_run')
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
        config={
            "delivery": "True",
            "rules_to_deliver": ",".join(rules_to_deliver)
        },
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
        delivery_json=delivery_json, sample_config_dict=sample_config_dict)

    # Add Housekeeper file to report
    delivery_json["files"].append({
        "path":
        report_file_name,
        "step":
        "balsamic_delivery",
        "format":
        get_file_extension(report_file_name),
        "tag": ["balsamic-report"],
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
        get_file_extension(sample_config),
        "tag": ["balsamic-config"],
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
        get_file_extension(sample_config_dict["analysis"]["dag"]),
        "tag": ["balsamic-dag"],
        "id":
        case_name,
    })
    # Add balsamic_qc_report
    if balsamic_qc_report:
        delivery_json["files"].append({
            "path":
            balsamic_qc_report,
            "step":
            "balsamic_delivery",
            "format":
            get_file_extension(balsamic_qc_report),
            "tag": ["delivery_report"],
            "id":
            case_name,
        })

    write_json(delivery_json, delivery_file_name)
    with open(delivery_file_name + ".yaml", "w") as fn:
        yaml.dump(delivery_json, fn, default_flow_style=False)

    LOG.info(f"Housekeeper delivery file {delivery_file_name}")
