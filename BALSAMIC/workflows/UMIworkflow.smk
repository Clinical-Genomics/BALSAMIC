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
table_dir = get_result_dir(config) + "/tables/"
plot_dir = get_result_dir(config) + "/plots/"
qc_dir = get_result_dir(config) + "/qc/"

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
VAR_CALLER = ["TNscope.umi","vardict.umi"]
ALL_STEPS = ["consensusaligned","consensusfiltered", "umialign"]
FILTERED_STEPS = ["consensusaligned","consensusfiltered"]
NEW_CASE_NAME = expand("{case_nm}.{step}", case_nm = CASE_NAME, step=FILTERED_STEPS)

# Define outputs
analysis_output = [ expand(vcf_dir + "SNV.somatic.{case_name}.{step}.{var_caller}.vcf.gz", case_name=CASE_NAME, step =FILTERED_STEPS, var_caller=VAR_CALLER), 
expand(vep_dir + "{var_type}.somatic.{case_name}.{var_caller}.{filters}.vcf.gz", var_type= "SNV", case_name=NEW_CASE_NAME, var_caller= VAR_CALLER, filters=["all", "pass"]),
expand(qc_dir + "{case_name}.{step}.umimetrics", case_name=CASE_NAME, step=FILTERED_STEPS),
expand(qc_dir + "{case_name}.{step}.collect_hsmetric_umi", case_name=CASE_NAME, step=FILTERED_STEPS),
expand(qc_dir + "{case_name}.{step}.mean_family_depth", case_name=CASE_NAME, step = FILTERED_STEPS),
expand(qc_dir + "{case_name}.{var_caller}.noiseAF", case_name=CASE_NAME, var_caller=['TNscope.umi']),
expand(plot_dir + "{case_name}.{var_caller}.AFplot.pdf", case_name=CASE_NAME, var_caller=['TNscope.umi']),
expand(table_dir + "{case_name}.{step}.{varcaller}.AFtable.txt", case_name=CASE_NAME, varcaller=VAR_CALLER, step= FILTERED_STEPS) ] 

config["rules"] = umi_call + variant_call + generate_tables + annotate_vcf + qc

for r in config["rules"]:
    include: Path(RULE_DIRECTORY, r).as_posix()

rule all:
    input:
        analysis_output
    output:
        os.path.join(get_result_dir(config), "analysis_finish")
    shell:
        "date +'%Y-%m-%d T%T %:z' > {output}"
