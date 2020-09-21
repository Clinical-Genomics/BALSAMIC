# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import os
import logging

from yapf.yapflib.yapf_api import FormatFile
from BALSAMIC.utils.cli import write_json
from BALSAMIC.utils.rule import get_rule_output 
from BALSAMIC.utils.rule import get_result_dir

LOG = logging.getLogger(__name__)

shell.prefix("set -eo pipefail; ")

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

singularity_image = config['singularity']['image'] 

# explictly check if cluster_config dict has zero keys.
if len(cluster_config.keys()) == 0:
    cluster_config = config

pre_align = ["snakemake_rules/quality_control/fastp.rule",
             "snakemake_rules/quality_control/fastqc.rule"]

align_qc = ["snakemake_rules/align/bwa_mem.rule",
            "snakemake_rules/quality_control/picard.rule",
            "snakemake_rules/quality_control/sambamba_depth.rule",
            "snakemake_rules/quality_control/mosdepth.rule",
            "snakemake_rules/quality_control/multiqc.rule",
            "snakemake_rules/quality_control/GATK.rule"]

config["rules"] = pre_align + align_qc

for r in config["rules"]:
    include: os.path.join(rule_dir + r)

if 'delivery' in config:
    wildcard_dict = { "sample": list(config["samples"].keys()),
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
            LOG.warning("Cannot deliver step (rule) {}: {}".format(my_rule,e))
            continue

        LOG.info("Delivering step (rule) {}.".format(my_rule))
        output_files_ready.extend(get_rule_output(rules=rules, rule_name=my_rule, output_file_wildcards=wildcard_dict))

    output_files_ready = [dict(zip(output_files_ready[0], value)) for value in output_files_ready[1:]]
    delivery_ready = os.path.join(get_result_dir(config),
                                  "delivery_report",
                                  config["analysis"]["case_id"] + "_delivery_ready.hk" )
    write_json(output_files_ready, delivery_ready)
    FormatFile(delivery_ready) 


rule all:
  input:
    os.path.join(*([result_dir + "qc/" + "multiqc_report.html"])),
  output:
    os.path.join(get_result_dir(config), "analysis_finish")
  shell:
    "date +'%Y-%m-%d T%T %:z' > {output}"

