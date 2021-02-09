# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import os
import logging

from BALSAMIC.utils.rule import (get_threads, get_result_dir,
                                 get_sample_type, get_script_path, get_vcf)
from BALSAMIC.utils.models import UMIworkflowConfig
from BALSAMIC.utils.constants import RULE_DIRECTORY, VCFANNO_TOML, umiworkflow_params
from BALSAMIC.utils.workflowscripts import get_densityplot

LOG = logging.getLogger(__name__)

shell.prefix("set -eo pipefail; ")

fastq_dir = get_result_dir(config) + "/fastq/"
benchmark_dir = config["analysis"]["benchmark"]
umi_dir = get_result_dir(config) + "/umi/"
vcf_dir = get_result_dir(config) + "/vcf/"
vep_dir = get_result_dir(config) + "/vep/"
umi_qc_dir = get_result_dir(config) + "/qc/"

singularity_image = config["singularity"]["image"]

# Declare sentieon variables
sentieon = True
SENTIEON_LICENSE = ''
SENTIEON_INSTALL_DIR = ''


# explictly check if cluster_config dict has zero keys.
if len(cluster_config.keys()) == 0:
    cluster_config = config

try:
    config["SENTIEON_LICENSE"] = os.environ["SENTIEON_LICENSE"]
    config["SENTIEON_EXEC"] = Path(os.environ["SENTIEON_INSTALL_DIR"], "bin", "sentieon").as_posix()
except Exception as error:
    LOG.error("ERROR: Set SENTIEON_LICENSE and SENTIEON_INSTALL_DIR environment variable to run this pipeline.")
    raise

# Define umiworkflow rules
umi_call = [
    "snakemake_rules/umi/sentieon_umiextract.rule",
    "snakemake_rules/umi/sentieon_consensuscall.rule"
]

variant_call = [
    "snakemake_rules/umi/sentieon_varcall_tnscope.rule",
    "snakemake_rules/umi/varcall_vardict.rule"
]

annotate_vcf = ["snakemake_rules/annotation/vep.rule"]

qc = ["snakemake_rules/umi/qc_umi.rule"]

generate_tables = ["snakemake_rules/umi/generate_AF_tables.rule"]

# parse parameters as workflow constants
paramsumi = UMIworkflowConfig.parse_obj(umiworkflow_params)

# Define wildcards
SAMPLES = config["samples"]
CASE_NAME = config["analysis"]["case_id"]
VAR_CALLER = ["TNscope","vardict"]
ALL_STEPS = ["consensusaligned","consensusfiltered", "umialign"]
FILTERED_STEPS = ["consensusaligned","consensusfiltered"]
EXTN_NAME = expand("{var_caller}_{step}_umi", var_caller = ["TNscope","vardict"], step=FILTERED_STEPS)

# Define outputs
analysis_output = [ expand(vep_dir + "{var_type}.somatic.{case_name}.{var_caller}.{filters}.vcf.gz", var_type= "SNV", case_name=CASE_NAME, var_caller=EXTN_NAME, filters=["all", "pass"]),
expand(umi_qc_dir + "{case_name}.{step}_umi.{metric}", case_name=CASE_NAME, step=FILTERED_STEPS, metric = ["metrics","collect_hsmetric", "mean_family_depth"]),
expand(umi_qc_dir + "{case_name}.{var_caller}_umi.noiseAF", case_name=CASE_NAME, var_caller=["TNscope"]),
expand(umi_qc_dir + "{case_name}.{var_caller}_umi.AFplot.pdf", case_name=CASE_NAME, var_caller=["TNscope"]),
expand(umi_qc_dir + "{case_name}.{var_caller}.AFtable.txt", case_name=CASE_NAME, var_caller=EXTN_NAME)
]

config["rules"] = umi_call + variant_call +  annotate_vcf + qc + generate_tables

for r in config["rules"]:
    include: Path(RULE_DIRECTORY, r).as_posix()

rule all:
    input:
        analysis_output
    output:
        os.path.join(get_result_dir(config), "analysis_finish")
    shell:
        "date +'%Y-%m-%d T%T %:z' > {output}"
