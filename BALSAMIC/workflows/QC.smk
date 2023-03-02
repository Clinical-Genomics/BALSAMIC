# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
import os
import logging
import tempfile

from pathlib import Path
from yapf.yapflib.yapf_api import FormatFile

from snakemake.exceptions import RuleException, WorkflowError

from BALSAMIC.utils.exc import BalsamicError

from BALSAMIC.utils.cli import (check_executable, generate_h5)
from BALSAMIC.utils.io import write_json

from BALSAMIC.utils.models import BalsamicWorkflowConfig

from BALSAMIC.utils.rule import (validate_fastq_input, get_fastqpatterns, get_mapping_info, get_fastq_info, get_rule_output, get_result_dir,
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
Path.mkdir(Path(tmp_dir), parents=True, exist_ok=True)

case_id = config["analysis"]["case_id"]
analysis_dir = config["analysis"]["analysis_dir"] + "/" + case_id + "/"
benchmark_dir = config["analysis"]["benchmark"]
fastqinput_dir = config["analysis"]["fastq_path"]
fastq_dir = get_result_dir(config) + "/fastq/"
bam_dir = get_result_dir(config) + "/bam/"
fastqc_dir = get_result_dir(config) + "/fastqc/"
result_dir = get_result_dir(config) + "/"
qc_dir = get_result_dir(config) + "/qc/"
delivery_dir = get_result_dir(config) + "/delivery/"

singularity_image = config['singularity']['image']


# Prepare sample_dict
sample_dict = {}
for sample in config["samples"]:
    sample_type = config["samples"][sample]["type"]
    if sample_type == "tumor":
        tumor_sample = sample
        sample_dict[tumor_sample] = get_fastq_info(tumor_sample, fastqinput_dir)
        sample_dict[tumor_sample]["sample_type"] = "TUMOR"
    else:
        normal_sample = sample
        sample_dict[normal_sample] = get_fastq_info(normal_sample, fastqinput_dir)
        sample_dict[normal_sample]["sample_type"] = "NORMAL"

# Validate fastq-info
validate_fastq_input(sample_dict, fastqinput_dir)

# Get fastq pattern --> fastq mapping
fastq_dict = {}
for sample in sample_dict:
    for fastq_pattern in sample_dict[sample]["fastqpair_patterns"]:
        fastq_dict[fastq_pattern] = sample_dict[sample]["fastqpair_patterns"][fastq_pattern]

# picarddup flag
picarddup = get_picard_mrkdup(config)

# Get mapping info
for sample in sample_dict:
    sample_dict[sample]["bam"] = get_mapping_info(samplename=sample,
                                    sample_dict=sample_dict,
                                    bam_dir=bam_dir,
                                    sequencing_type=config["analysis"]["sequencing_type"])


# parse parameters as constants to workflows
params = BalsamicWorkflowConfig.parse_obj(WORKFLOW_PARAMS)

# Capture kit name
if config["analysis"]["sequencing_type"] != "wgs":
    capture_kit = os.path.split(config["panel"]["capture_kit"])[1]

# Set case id/name
case_id = config["analysis"]["case_id"]

# explicitly check if cluster_config dict has zero keys.
if len(cluster_config.keys()) == 0:
    cluster_config = config

if "hg38" in config["reference"]["reference_genome"]:
    config["reference"]["genome_version"] = "hg38"
elif "canfam3" in config["reference"]["reference_genome"]:
    config["reference"]["genome_version"] = "canfam3"
else:
    config["reference"]["genome_version"] = "hg19"

LOG.info('Genome version set to %s', config["reference"]["genome_version"])


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
                "snakemake_rules/quality_control/qc_metrics.rule"
]

if "paired" in config['analysis']['analysis_type']:
    rules_to_include.append("snakemake_rules/variant_calling/mergetype_normal.rule")

    # Somalier only implemented for hg38 and hg19
    if "canfam3" not in config["reference"]["reference_genome"]:
        rules_to_include.append("snakemake_rules/quality_control/somalier.rule")

for r in rules_to_include:
    include: Path(RULE_DIRECTORY, r).as_posix()
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

        # Delete a temporal directory tree
        try:
            shutil.rmtree(params.tmp_dir)
        except OSError as e:
            print ("Error: %s - %s." % (e.filename, e.strerror))

        # Finish timestamp file
        with open(str(output.finish_file), mode="w") as finish_file:
            finish_file.write("%s\n" % datetime.datetime.now())
