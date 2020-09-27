# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import os
import logging

from yapf.yapflib.yapf_api import FormatFile

from snakemake.exceptions import RuleException, WorkflowError
from BALSAMIC.utils.cli import write_json
from BALSAMIC.utils.rule import get_variant_callers
from BALSAMIC.utils.rule import get_rule_output
from BALSAMIC.utils.rule import get_result_dir
from BALSAMIC.utils.rule import get_vcf

shell.prefix("set -eo pipefail; ")

LOG = logging.getLogger(__name__)

# Set temporary dir environment variable
os.environ['TMPDIR'] = get_result_dir(config)

tmp_dir = os.path.join(get_result_dir(config), "tmp")
rule_dir = config["rule_directory"]
benchmark_dir = config["analysis"]["benchmark"]
fastq_dir = get_result_dir(config) + "/fastq/"
bam_dir = get_result_dir(config) + "/bam/"
cnv_dir = get_result_dir(config) + "/cnv/"
cutadapt_dir = get_result_dir(config) + "/cutadapt/"
fastqc_dir = get_result_dir(config) + "/fastqc/"
result_dir = get_result_dir(config) + "/"
vcf_dir = get_result_dir(config) + "/vcf/"
vep_dir = get_result_dir(config) + "/vep/"
qc_dir = result_dir + "qc/"
delivery_dir = get_result_dir(config) + "/delivery/"

singularity_image = config['singularity']['image']

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
    sentieon = False
    LOG.warn("Set environment variables SENTIEON_LICENSE and SENTIEON_INSTALL_DIR to run SENTIEON variant callers")

# Define set of rules
qc_rules = [
  "snakemake_rules/quality_control/fastp.rule",
  "snakemake_rules/quality_control/fastqc.rule",
  "snakemake_rules/quality_control/GATK.rule",
  "snakemake_rules/quality_control/multiqc.rule",
  "snakemake_rules/quality_control/picard.rule",
  "snakemake_rules/quality_control/sambamba_depth.rule",
  "snakemake_rules/quality_control/mosdepth.rule"
  ]

align_rules = [
  "snakemake_rules/align/bwa_mem.rule"
  ]

annotation_rules = [
  "snakemake_rules/annotation/vep.rule"
  ]

variantcalling_rules = [
  "snakemake_rules/variant_calling/germline.rule",
  "snakemake_rules/variant_calling/split_bed.rule"
  ]

germline_caller = ["haplotypecaller", "strelka_germline", "manta_germline"]

if sentieon:
    germline_caller.append('dnascope')

if config['analysis']['analysis_type'] == "paired":

    qc_rules.append("snakemake_rules/quality_control/contest.rule")

    variantcalling_rules.extend([
      "snakemake_rules/variant_calling/somatic_tumor_normal.rule",
      "snakemake_rules/variant_calling/somatic_sv_tumor_normal.rule",
      "snakemake_rules/variant_calling/mergetype.rule",
      "snakemake_rules/variant_calling/cnvkit_paired.rule"
      ])

    somatic_caller_snv = get_variant_callers(config=config, analysis_type="paired", workflow_solution="BALSAMIC", mutation_type="SNV", mutation_class="somatic")
    sentieon_callers = ["tnhaplotyper"] if sentieon else []
    somatic_caller_sv = ["manta", "cnvkit"]

else:

    annotation_rules.append("snakemake_rules/annotation/varcaller_filter.rule")

    variantcalling_rules.extend([
      "snakemake_rules/variant_calling/cnvkit_single.rule",
      "snakemake_rules/variant_calling/mergetype_tumor.rule",
      "snakemake_rules/variant_calling/somatic_tumor_only.rule",
      "snakemake_rules/variant_calling/somatic_sv_tumor_only.rule"
      ])

    somatic_caller_snv = get_variant_callers(config=config, analysis_type="single", workflow_solution="BALSAMIC", mutation_type="SNV", mutation_class="somatic")
    sentieon_callers = ["tnhaplotyper"] if sentieon else []
    somatic_caller_sv = ["manta", "cnvkit"]

somatic_caller = somatic_caller_snv + somatic_caller_sv + sentieon_callers
if "disable_variant_caller" in config:
    somatic_caller.remove(config["disable_variant_caller"])

config["rules"] = align_rules + qc_rules


for r in config["rules"]:
    include: os.path.join(rule_dir + r)

# Define common and analysis specific outputs
quality_control_results = [ result_dir + "qc/" + "multiqc_report.html" ]

analysis_specific_results = []
if config['analysis']["analysis_type"] in ["paried", "single"]:
    config["rules"] = config["rules"] + variantcalling_rules + annotation_rules
    analysis_specific_results = [expand(vep_dir + "{vcf}.vcf.gz",
                                        vcf=get_vcf(config, germline_caller, config["samples"])),
                                 expand(vep_dir + "{vcf}.{filters}.vcf.gz",
                                        vcf=get_vcf(config, somatic_caller, [config["analysis"]["case_id"]]),
                                        filters = ["all", "pass"]),
                                 expand(vep_dir + "{vcf}.pass.balsamic_stat",
                                        vcf=get_vcf(config, ["vardict"], [config["analysis"]["case_id"]]))]

if config['analysis']['analysis_type'] == "single":
    analysis_specific_results.extend(expand(vep_dir + "{vcf}.all.filtered.vcf.gz",
                                            vcf=get_vcf(config, ["vardict"], [config["analysis"]["case_id"]])))

if 'delivery' in config:
    wildcard_dict = { "sample": list(config["samples"].keys()),
                      "case_name": config["analysis"]["case_id"],
                      "allow_missing": True
                    }

    if config['analysis']["analysis_type"] in ["paried", "single"]:
        wildcard_dict.update({"var_type": ["CNV", "SNV", "SV"],
                              "var_class": ["somatic", "germline"],
                              "var_caller": somatic_caller + germline_caller,
                              "bedchrom": config["panel"]["chrom"] if "panel" in config else [],
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
            LOG.warning("Cannot deliver step (rule) {}: {}".format(my_rule,e))
            continue

        LOG.info("Delivering step (rule) {}.".format(my_rule))
        output_files_ready.extend(get_rule_output(rules=rules,
                                                  rule_name=my_rule,
                                                  output_file_wildcards=wildcard_dict))

    output_files_ready = [dict(zip(output_files_ready[0], value)) for value in output_files_ready[1:]]
    delivery_ready = os.path.join(get_result_dir(config),
                                  "delivery_report",
                                  config["analysis"]["case_id"] + "_delivery_ready.hk" )
    write_json(output_files_ready, delivery_ready)
    FormatFile(delivery_ready)

rule all:
    input:
        quality_control_results,
        analysis_specific_results
    output:
        os.path.join(get_result_dir(config), "analysis_finish")
    shell:
        "date +'%Y-%m-%d T%T %:z' > {output}"
