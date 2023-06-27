# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from pathlib import Path
import glob
import tempfile
import os
import logging


from BALSAMIC.utils.exc import BalsamicError

from BALSAMIC.utils.rule import get_fastqpatterns, get_mapping_info, get_picard_mrkdup, get_threads, get_result_dir
from BALSAMIC.constants.common import RULE_DIRECTORY
from BALSAMIC.constants.workflow_params import WORKFLOW_PARAMS
from BALSAMIC.utils.models import BalsamicWorkflowConfig

shell.prefix("set -eo pipefail; ")

localrules: all

LOG = logging.getLogger(__name__)
logging.getLogger("filelock").setLevel("WARN")

# parse parameters as constants to workflows
params = BalsamicWorkflowConfig.parse_obj(WORKFLOW_PARAMS)

fastq_dir = get_result_dir(config) + "/fastq/"
analysis_dir = get_result_dir(config)
qc_dir = analysis_dir + "/qc/"
bam_dir =  analysis_dir + "/bam/"
cnv_dir =  analysis_dir + "/cnv/"

reffasta = config["reference"]["reference_genome"]
refflat = config["reference"]["refflat"]
access_5kb_hg19 = config["reference"]["access_regions"]
target_bed = config["panel"]["capture_kit"]
singularity_image = config["singularity"]["image"]
benchmark_dir = config["analysis"]["benchmark"]
version = config["analysis"]["pon_version"]

tmp_dir = os.path.join(analysis_dir, "tmp", "" )
Path.mkdir(Path(tmp_dir), parents=True, exist_ok=True)

# Prepare sample_dict
sample_dict = dict(config["samples"])
for sample in sample_dict:
    sample_type = sample_dict[sample]["type"]
    if sample_type == "tumor":
        tumor_sample = sample
        sample_dict[tumor_sample]["sample_type"] = "TUMOR"
    else:
        normal_sample = sample
        sample_dict[normal_sample]["sample_type"] = "NORMAL"


# Get fastq pattern --> fastq mapping
fastq_dict = {}
for sample in sample_dict:
    for fastq_pattern in sample_dict[sample]["fastq_info"]:
        fastq_dict[fastq_pattern] = sample_dict[sample]["fastq_info"][fastq_pattern]

# Get mapping info
for sample in sample_dict:
    sample_dict[sample]["bam"] = get_mapping_info(samplename=sample,
                                    sample_dict=sample_dict,
                                    bam_dir=bam_dir,
                                    analysis_type=config["analysis"]["analysis_type"])

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

if not Path(config["SENTIEON_EXEC"]).exists():
    LOG.error("Sentieon executable not found {}".format(Path(config["SENTIEON_EXEC"]).as_posix()))
    raise BalsamicError

panel_name = os.path.split(target_bed)[1].replace('.bed','')

coverage_references = expand(cnv_dir + "{sample}.{cov}coverage.cnn", sample=config["samples"], cov=['target','antitarget'])
baited_beds = expand(cnv_dir + "{cov}.bed", cov=['target','antitarget'])
pon_reference = expand(cnv_dir + panel_name + "_CNVkit_PON_reference_" + version + ".cnn")
pon_finish = expand(analysis_dir + "/" + "analysis_PON_finish")

config["rules"] = [
    "snakemake_rules/quality_control/fastp.rule",
    "snakemake_rules/align/sentieon_alignment.rule",
]

for r in config["rules"]:
    include: Path(RULE_DIRECTORY, r).as_posix()

rule all:
    input:
        ref_cnn = pon_reference
    output:
        pon_finish_file = pon_finish
    run:
        import datetime

        # PON finish timestamp file
        with open(str(output.pon_finish_file), mode="w") as finish_file:
            finish_file.write("%s\n" % datetime.datetime.now())

rule create_target:
    input:
        target_bait = target_bed,
        refFlat = refflat,
        access_bed = access_5kb_hg19
    output:
        target_bed = cnv_dir + "target.bed",
        offtarget_bed = cnv_dir + "antitarget.bed"
    singularity:
        Path(singularity_image, "varcall_cnvkit.sif").as_posix()
    benchmark:
        Path(benchmark_dir, "cnvkit.targets.tsv").as_posix()
    shell:
        """
cnvkit.py target {input.target_bait} --annotate {input.refFlat} --split -o {output.target_bed};
cnvkit.py antitarget {input.target_bait} -g {input.access_bed} -o {output.offtarget_bed};
        """

rule create_coverage:
    input:
        bam = lambda wildcards: sample_dict[wildcards.sample]["bam"]["final_bam"],
        target_bed = cnv_dir + "target.bed",
        antitarget_bed = cnv_dir + "antitarget.bed"
    output:
        target_cnn = cnv_dir + "{sample}.targetcoverage.cnn",
        antitarget_cnn = cnv_dir + "{sample}.antitargetcoverage.cnn"
    singularity:
        Path(singularity_image, "varcall_cnvkit.sif").as_posix()
    benchmark:
        Path(benchmark_dir, "cnvkit_{sample}.coverage.tsv").as_posix()
    shell:
        """
cnvkit.py coverage {input.bam} {input.target_bed} -o {output.target_cnn};
cnvkit.py coverage {input.bam} {input.antitarget_bed} -o {output.antitarget_cnn};
        """

rule create_reference:
    input:
        cnn = expand(cnv_dir + "{sample}.{prefix}coverage.cnn", sample=config["samples"], prefix=["target", "antitarget"]),
        ref = reffasta
    output:
        ref_cnn = pon_reference
    singularity:
        Path(singularity_image, "varcall_cnvkit.sif").as_posix()
    benchmark:
        Path(benchmark_dir, "cnvkit.reference.tsv").as_posix()
    shell:
        """
cnvkit.py reference {input.cnn} --fasta {input.ref} -o {output.ref_cnn} ;
        """
