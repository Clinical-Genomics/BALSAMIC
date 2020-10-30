# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import os
import logging

from BALSAMIC.utils.rule import get_result_dir
from BALSAMIC.utils.constants import RULE_DIRECTORY

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
    config["SENTIEON_INSTALL_DIR"] = os.environ["SENTIEON_INSTALL_DIR"]
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

annotate_vcf = ["snakemake_rules/umi/annotate_vep.rule"]

qc = ["snakemake_rules/umi/qc_umi.rule"]

generate_tables = ["snakemake_rules/umi/generate_AF_tables.rule"]

# Define wildcards
SAMPLES = config["samples"]
VAR_CALLER = ["TNscope","vardict"]
STEPS = ["umialign","consensusfiltered"]
#STEPS = ["consensusfiltered"]

# Define outputs
analysis_output = [ expand(vcf_dir + "{sample}.{var_caller}.{step}.vcf.gz", sample=SAMPLES, var_caller=VAR_CALLER, step = STEPS),
expand(vep_dir + "{sample}.{var_caller}.{step}.{filler}.vcf.gz", sample=SAMPLES, var_caller=VAR_CALLER, filler=["all","pass"], step=STEPS),
expand(qc_dir + "{sample}.{step}.umimetrics", sample=SAMPLES, step=STEPS),
expand(qc_dir + "{sample}.{step}.collect_hsmetric_umi", sample=SAMPLES, step=STEPS),
expand(qc_dir + "{sample}.{step}.mean_family_depth", sample=SAMPLES, step = STEPS),
expand(qc_dir + "{sample}.TNscope.noiseAF", sample=SAMPLES),
#expand(plot_dir + "{sample}.TNscope.AFplot.pdf", sample=SAMPLES),
expand(table_dir + "{sample}.{varcaller}.consensusfiltered.AFtable.txt", sample=SAMPLES, varcaller=VAR_CALLER) ]

config["rules"] = umi_call + variant_call + annotate_vcf + generate_tables + qc

for r in config["rules"]:
    include: Path(RULE_DIRECTORY, r).as_posix()

rule all:
    input:
        analysis_output
    output:
        os.path.join(get_result_dir(config), "analysis_finish")
    shell:
        "date +'%Y-%m-%d T%T %:z' > {output}"
