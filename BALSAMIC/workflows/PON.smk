# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from pathlib import Path
import glob
import tempfile
import os

from BALSAMIC.constants.paths import BALSAMIC_DIR
from BALSAMIC.utils.rule import get_picard_mrkdup, get_threads, get_result_dir
from BALSAMIC.constants.workflow_params import WORKFLOW_PARAMS
from BALSAMIC.models.models import BalsamicWorkflowConfig

shell.prefix("set -eo pipefail; ")

localrules: all

# parse parameters as constants to workflows
params = BalsamicWorkflowConfig.parse_obj(WORKFLOW_PARAMS)

fastq_dir =  config["analysis"]["fastq_path"]
analysis_dir = get_result_dir(config)
analysis_fastq_dir = analysis_dir + "/fastq/"
concat_dir = get_result_dir(config) + "/concat/"
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

picarddup = get_picard_mrkdup(config)

panel_name = os.path.split(target_bed)[1].replace('.bed','')

coverage_references = expand(cnv_dir + "{sample}.{cov}coverage.cnn", sample=config["samples"], cov=['target','antitarget'])
baited_beds = expand(cnv_dir + "{cov}.bed", cov=['target','antitarget'])
pon_reference = expand(cnv_dir + panel_name + "_CNVkit_PON_reference_" + version + ".cnn")
pon_finish = expand(analysis_dir + "/" + "analysis_PON_finish")

config["rules"] = [
    "snakemake_rules/concatenation/concatenation.rule",
    "snakemake_rules/quality_control/fastp.rule",
    "snakemake_rules/align/bwa_mem.rule",
]

for r in config["rules"]:
    include: Path(BALSAMIC_DIR, r).as_posix()

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
        bam = bam_dir + "{sample}.sorted." + picarddup  + ".bam",
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
