# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import os
import logging

from BALSAMIC.utils.rule import (get_threads, get_result_dir,
                                 get_sample_type, get_script_path, get_vcf)
from BALSAMIC.utils.models import UMIworkflowConfig
from BALSAMIC.utils.constants import RULE_DIRECTORY, VCFANNO_TOML, umiworkflow_params

LOG = logging.getLogger(__name__)

shell.prefix("set -eo pipefail; ")

fastq_dir = get_result_dir(config) + "/fastq/"
benchmark_dir = config["analysis"]["benchmark"]
umi_dir = get_result_dir(config) + "/umi/"
vcf_dir = get_result_dir(config) + "/vcf/"
vep_dir = get_result_dir(config) + "/vep/"
umi_qc_dir = get_result_dir(config) + "/qc/umi_qc/"
qc_dir = get_result_dir(config) + "/qc/umi_qc/"
tmp_dir = os.path.join(get_result_dir(config), "tmp", "" )
Path.mkdir(Path(tmp_dir), exist_ok=True)

singularity_image = config["singularity"]["image"]

# Declare sentieon variables
sentieon = True
SENTIEON_LICENSE = ''
SENTIEON_INSTALL_DIR = ''


# explictly check if cluster_config dict has zero keys.
if len(cluster_config.keys()) == 0:
    cluster_config = config

try:
    config["SENTIEON_INSTALL_DIR"] = os.environ["SENTIEON_INSTALL_DIR"]
    config["SENTIEON_LICENSE"] = os.environ["SENTIEON_LICENSE"]
    config["SENTIEON_EXEC"] = Path(os.environ["SENTIEON_INSTALL_DIR"], "bin", "sentieon").as_posix()
except Exception as error:
    LOG.error("ERROR: Set SENTIEON_LICENSE and SENTIEON_INSTALL_DIR environment variable to run this pipeline.")
    raise

# Sample names for tumor or normal
tumor_sample = get_sample_type(config["samples"], "tumor")[0]
if config['analysis']['analysis_type'] == "paired":
     normal_sample = get_sample_type(config["samples"], "normal")[0]

# Define umiworkflow rules
fastp_umi = ["snakemake_rules/quality_control/fastp.rule"]

umi_call = [
    "snakemake_rules/umi/sentieon_umiextract.rule",
    "snakemake_rules/umi/sentieon_consensuscall.rule",
    "snakemake_rules/umi/mergetype_tumor_umi.rule"
]

if config["analysis"]["analysis_type"] == "single":
    variant_call = ["snakemake_rules/umi/sentieon_varcall_tnscope.rule"]
else:
    variant_call = ["snakemake_rules/umi/sentieon_varcall_tnscope_tn.rule"]
    umi_call.extend(["snakemake_rules/umi/mergetype_normal_umi.rule"])

annotate_vcf = ["snakemake_rules/annotation/vep.rule"]

qc = ["snakemake_rules/umi/qc_umi.rule"]

generate_tables = ["snakemake_rules/umi/generate_AF_tables.rule"]

# parse parameters as workflow constants
paramsumi = UMIworkflowConfig.parse_obj(umiworkflow_params)

# Define wildcards
SAMPLES = config["samples"]
CASE_NAME = config["analysis"]["case_id"]

# Define outputs
analysis_output = [expand(vep_dir + "{var_type}.somatic.{case_name}.{var_caller}.all.vcf.gz", var_type= "SNV", case_name=CASE_NAME, var_caller=["TNscope_umi"]),
expand(umi_qc_dir + "{sample}.umi.{metric}", sample=SAMPLES, metric = ["metrics", "mean_family_depth"])]

config["rules"] = fastp_umi + umi_call + variant_call +  annotate_vcf + qc

if "background_variants" in config:
    analysis_output.extend([expand(umi_qc_dir + "{case_name}.{var_caller}.AFtable.txt",
                                   case_name = config["analysis"]["case_id"],
                                   var_caller =["TNscope_umi"])])
    config["rules"] = config["rules"] + generate_tables

for r in config["rules"]:
    include: Path(RULE_DIRECTORY, r).as_posix()

rule all:
    input:
        analysis_output
    output:
        os.path.join(get_result_dir(config), "analysis_finish")
    shell:
        "date +'%Y-%m-%d T%T %:z' > {output}"
