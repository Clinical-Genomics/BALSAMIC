# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
import os
import logging
import tempfile

from pathlib import Path
from typing import List
from yapf.yapflib.yapf_api import FormatFile

from snakemake.exceptions import RuleException, WorkflowError

from BALSAMIC.constants.rules import SNAKEMAKE_RULES
from BALSAMIC.constants.paths import BALSAMIC_DIR
from BALSAMIC.constants.analysis import FastqName, SampleType
from BALSAMIC.utils.exc import BalsamicError

from BALSAMIC.utils.cli import check_executable, generate_h5
from BALSAMIC.utils.io import write_json, write_finish_file

from BALSAMIC.models.analysis import BalsamicWorkflowConfig, BalsamicConfigModel

from BALSAMIC.utils.rule import (get_rule_output, get_result_dir,
                                 get_script_path, get_threads,
                                 get_sequencing_type, get_capture_kit)

from BALSAMIC.constants.workflow_params import WORKFLOW_PARAMS

# Initialize BalsamicConfigModel
config_model = BalsamicConfigModel.parse_obj(config)

shell.executable("/bin/bash")
shell.prefix("set -eo pipefail; ")

LOG = logging.getLogger(__name__)

# Get case id/name
case_id: str = config_model.analysis.case_id
# Get analysis dir
analysis_dir_home: str = config_model.analysis.analysis_dir
analysis_dir: str = os.path.join(analysis_dir_home, "analysis", case_id, "")
# Get result dir
result_dir: str = os.path.join(config_model.analysis.result, "")

# Create a temporary directory with trailing /
tmp_dir: str = os.path.join(result_dir, "tmp", "")
Path.mkdir(Path(tmp_dir), parents=True, exist_ok=True)

# Directories
input_fastq_dir: str = config_model.analysis.fastq_path + "/"
benchmark_dir: str = config_model.analysis.benchmark
fastq_dir: str = os.path.join(result_dir, "fastq", "")
bam_dir: str = os.path.join(result_dir, "bam", "")
fastqc_dir: str = os.path.join(result_dir, "fastqc", "")
vcf_dir: str = os.path.join(result_dir, "vcf", "")
qc_dir: str = os.path.join(result_dir, "qc", "")
delivery_dir: str = os.path.join(result_dir, "delivery", "")


# Run information
singularity_image: str = config_model.singularity['image']
sample_names: List[str] = config_model.get_all_sample_names()
tumor_sample: str = config_model.get_sample_name_by_type(SampleType.TUMOR)
if config_model.analysis.analysis_type == "paired":
    normal_sample: str = config_model.get_sample_name_by_type(SampleType.NORMAL)

# parse parameters as constants to workflows
params = BalsamicWorkflowConfig.parse_obj(WORKFLOW_PARAMS)

# Capture kit name
if config["analysis"]["sequencing_type"] != "wgs":
    capture_kit = os.path.split(config["panel"]["capture_kit"])[1]

# explicitly check if cluster_config dict has zero keys.
if len(cluster_config.keys()) == 0:
    cluster_config = config

# Find and set Sentieon binary and license server from env variables
try:
    config["SENTIEON_LICENSE"] = os.environ["SENTIEON_LICENSE"]
    config["SENTIEON_INSTALL_DIR"] = os.environ["SENTIEON_INSTALL_DIR"]

    if os.getenv("SENTIEON_EXEC") is not None:
        config["SENTIEON_EXEC"] = os.environ["SENTIEON_EXEC"]
    else:
        config["SENTIEON_EXEC"] = Path(os.environ["SENTIEON_INSTALL_DIR"], "bin", "sentieon").as_posix()
except KeyError as error:
    LOG.error("Set environment variables SENTIEON_LICENSE, SENTIEON_INSTALL_DIR, SENTIEON_EXEC "
              "to run SENTIEON variant callers")
    raise BalsamicError

if "hg38" in config["reference"]["reference_genome"]:
    config["reference"]["genome_version"] = "hg38"
elif "canfam3" in config["reference"]["reference_genome"]:
    config["reference"]["genome_version"] = "canfam3"
else:
    config["reference"]["genome_version"] = "hg19"

LOG.info('Genome version set to %s', config["reference"]["genome_version"])


# Set temporary dir environment variable
os.environ['TMPDIR'] = get_result_dir(config)

# Include rules
analysis_type = config['analysis']["analysis_type"]
sequence_type = config['analysis']["sequencing_type"]

rules_to_include = []
for workflow_type, value in SNAKEMAKE_RULES.items():
    if workflow_type in ["common", analysis_type + "_" + sequence_type]:
        rules_to_include.extend(value.get("qc", []) + value.get("align", []))
rules_to_include = [rule for rule in rules_to_include if "umi" not in rule]
if "snakemake_rules/quality_control/report.rule" in rules_to_include:
    rules_to_include = [rule for rule in rules_to_include if "quality_control/report.rule" not in rule]

# Somalier only implemented for hg38 and hg19
if "canfam3" in config["reference"]["reference_genome"]:
    rules_to_include.remove("snakemake_rules/quality_control/somalier.rule")

for r in rules_to_include:
    include: Path(BALSAMIC_DIR, r).as_posix()

LOG.info(f"The following rules will be included in the workflow: {rules_to_include}")

# Define common and analysis specific outputs
quality_control_results = [
    os.path.join(qc_dir, case_id + "_metrics_deliverables.yaml"),
    os.path.join(qc_dir, "multiqc_report.html"),
]

if 'delivery' in config:
    wildcard_dict = {"sample": list(config["samples"].keys())+["tumor", "normal"],
                     "case_name": config["analysis"]["case_id"],
                     "allow_missing": True
                     }

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
        finish_file = os.path.join(get_result_dir(config), "analysis_finish")
    params:
        tmp_dir = tmp_dir,
    run:
        import datetime
        import shutil

        # Remove temporary directory tree
        try:
            shutil.rmtree(params.tmp_dir)
        except OSError as e:
            print ("Error: %s - %s." % (e.filename, e.strerror))

        # Finish timestamp file
        write_finish_file(file_path=output.finish_file)
