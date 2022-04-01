# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
import os
import logging
import tempfile

from pathlib import Path
from yapf.yapflib.yapf_api import FormatFile

from snakemake.exceptions import RuleException, WorkflowError

from BALSAMIC.utils.exc import BalsamicError

from BALSAMIC.utils.cli import (write_json, check_executable, generate_h5)

from BALSAMIC.utils.models import BalsamicWorkflowConfig

from BALSAMIC.utils.rule import (get_rule_output, get_result_dir,
                                 get_sample_type, get_picard_mrkdup, get_script_path,
                                 get_threads, get_sequencing_type, get_capture_kit)

from BALSAMIC.constants.common import (RULE_DIRECTORY);
from BALSAMIC.constants.workflow_params import WORKFLOW_PARAMS


shell.executable("/bin/bash")
shell.prefix("set -eo pipefail; ")

LOG = logging.getLogger(__name__)
logging.getLogger("filelock").setLevel("WARN")

# Create a temporary directory with trailing /
tmp_dir = os.path.join(get_result_dir(config), "tmp", "" )
Path.mkdir(Path(tmp_dir), exist_ok=True)

benchmark_dir = config["analysis"]["benchmark"]
fastq_dir = get_result_dir(config) + "/fastq/"
bam_dir = get_result_dir(config) + "/bam/"
fastqc_dir = get_result_dir(config) + "/fastqc/"
result_dir = get_result_dir(config) + "/"
qc_dir = get_result_dir(config) + "/qc/"
delivery_dir = get_result_dir(config) + "/delivery/"

singularity_image = config['singularity']['image']

# picarddup flag
picarddup = get_picard_mrkdup(config)

# parse parameters as constants to workflows
params = BalsamicWorkflowConfig.parse_obj(WORKFLOW_PARAMS)

# Capture kit name
if config["analysis"]["sequencing_type"] != "wgs":
    capture_kit = os.path.split(config["panel"]["capture_kit"])[1]

# Sample names for tumor or normal
tumor_sample = get_sample_type(config["samples"], "tumor")[0]
if "paired" in config['analysis']['analysis_type']:
    normal_sample = get_sample_type(config["samples"], "normal")[0]

# Set case id/name
case_id = config["analysis"]["case_id"]

# explicitly check if cluster_config dict has zero keys.
if len(cluster_config.keys()) == 0:
    cluster_config = config

# Add reference assembly if not defined for backward compatibility
if 'genome_version' not in config["reference"]:
    GENOME_VERSION = 'hg19' ## if hg19 convention works, replace accordingly
    LOG.info('Genome version was not found in config. Setting it to %s', GENOME_VERSION)

# Add normal sample if analysis is paired
germline_call_samples = ["tumor"]
if config['panel']:
    germline_call_samples.append("normal")

# Create list of chromosomes in panel for panel only variant calling to be used in rules
if config["analysis"]["sequencing_type"] != "wgs":
    chromlist = config["panel"]["chrom"]

# Set temporary dir environment variable
os.environ['TMPDIR'] = get_result_dir(config)

analysis_type = config['analysis']["analysis_type"]

rules_to_include = [
                "snakemake_rules/quality_control/fastp.rule",
                "snakemake_rules/quality_control/fastqc.rule",
                "snakemake_rules/quality_control/multiqc.rule",
                "snakemake_rules/variant_calling/mergetype_tumor.rule",
                "snakemake_rules/quality_control/picard.rule",
                "snakemake_rules/quality_control/sambamba_depth.rule",
                "snakemake_rules/quality_control/mosdepth.rule",
                "snakemake_rules/align/bwa_mem.rule",
]
if "paired" in config['analysis']['analysis_type']:
    rules_to_include.append("snakemake_rules/variant_calling/mergetype_normal.rule")

# for r in rules_to_include:
for r in rules_to_include:
    include: Path(RULE_DIRECTORY, r).as_posix()
LOG.info(f"The following rules will be included in the workflow: {rules_to_include}")

# Define common and analysis specific outputs
quality_control_results = [result_dir + "qc/" + "multiqc_report.html"]

if 'delivery' in config:
    wildcard_dict = {"sample": list(config["samples"].keys())+["tumor", "normal"],
                     "case_name": config["analysis"]["case_id"],
                     "allow_missing": True
                     }

    if config['analysis']["analysis_type"] in ["paired", "single"]:
        wildcard_dict.update({"bedchrom": config["panel"]["chrom"] if "panel" in config else [],
                              })

    if 'rules_to_deliver' in config:
        rules_to_deliver = config['rules_to_deliver'].split(",")
    else:
        rules_to_deliver = ['multiqc']

    output_files_ready = [('path', 'path_index', 'step', 'tag', 'id', 'format')]

    for my_rule in set(rules_to_deliver):
        try:
            housekeeper_id = getattr(rules, my_rule).params.housekeeper_id
        except (ValueError, AttributeError, RuleException, WorkflowError) as e:
            LOG.warning("Cannot deliver step (rule) {}: {}".format(my_rule, e))
            continue

        LOG.info("Delivering step (rule) {} {}.".format(my_rule, housekeeper_id))
        files_to_deliver = get_rule_output(rules=rules, rule_name=my_rule, output_file_wildcards=wildcard_dict)
        LOG.debug("The following files added to delivery: {}".format(files_to_deliver))
        output_files_ready.extend(files_to_deliver)

    output_files_ready = [dict(zip(output_files_ready[0], value)) for value in output_files_ready[1:]]
    delivery_ready = os.path.join(get_result_dir(config),
                                  "delivery_report",
                                  config["analysis"]["case_id"] + "_delivery_ready.hk")
    write_json(output_files_ready, delivery_ready)
    FormatFile(delivery_ready)

rule all:
    input:
        quality_control_results
    output:
        qc_json_file = os.path.join(get_result_dir(config), "qc", "qc_metrics_summary.json"),
        finish_file = os.path.join(get_result_dir(config), "analysis_finish")
    params:
        tmp_dir = tmp_dir,
        result_dir = result_dir,
        sequencing_type = get_sequencing_type(config),
        panel_bed = get_capture_kit(config)
    run:
        import datetime
        import shutil

        from BALSAMIC.utils.qc_metrics import get_qc_metrics_json

        # Save QC metrics to a JSON file
        try:
            qc_metrics_summary = get_qc_metrics_json(params.result_dir, params.sequencing_type, params.panel_bed)
            write_json(qc_metrics_summary, str(output.qc_json_file))
        except ValueError as val_exc:
            LOG.error(val_exc)
            raise BalsamicError

        # Delete a temporal directory tree
        try:
            shutil.rmtree(params.tmp_dir)
        except OSError as e:
            print ("Error: %s - %s." % (e.filename, e.strerror))

        # Finish timestamp file
        with open(str(output.finish_file), mode="w") as finish_file:
            finish_file.write("%s\n" % datetime.datetime.now())
