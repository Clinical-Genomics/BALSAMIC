# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from pathlib import Path
import glob
import tempfile
import os


from BALSAMIC.utils.rule import (get_picard_mrkdup, get_threads, 
                                 get_result_dir, get_pon_samples)
from BALSAMIC.utils.constants import RULE_DIRECTORY

shell.prefix("set -eo pipefail; ")

localrules: all

fastq_dir = get_result_dir(config) + "/fastq/"
qc_dir = get_result_dir(config) + "/qc/"
bam_dir =  get_result_dir(config) + "/bam/"
cnv_dir =  get_result_dir(config) + "/cnv/"

reffasta = config["reference"]["reference_genome"]
refflat = config["reference"]["refflat"]
access_5kb_hg19 = config["reference"]["access_regions"]
target_bed = config["panel"]["capture_kit"]
singularity_image = config["singularity"]["image"]
benchmark_dir = config["analysis"]["benchmark"]

tmp_dir = os.path.join(get_result_dir(config), "tmp", "" )
Path.mkdir(Path(tmp_dir), exist_ok=True)

picarddup = get_picard_mrkdup(config)
samples = get_pon_samples(fastq_dir)


ALL_COVS = expand(cnv_dir + "{sample}.{cov}coverage.cnn", sample=samples, cov=['target','antitarget'])
ALL_REFS = expand(cnv_dir + "{cov}.bed", cov=['target','antitarget'])
ALL_PON = expand(cnv_dir + config["analysis"]["case_id"] + "_PON_reference.cnn")
PON_DONE = expand(cnv_dir + "PON." + "reference" + ".done")

config["rules"] = ["snakemake_rules/quality_control/fastp.rule", 
                   "snakemake_rules/align/bwa_mem.rule"]

for r in config["rules"]:
    include: Path(RULE_DIRECTORY, r).as_posix()

rule all:
    input: ALL_REFS + ALL_COVS +  PON_DONE

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
source activate varcall_cnvkit;
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
source activate varcall_cnvkit;
cnvkit.py coverage {input.bam} {input.target_bed} -o {output.target_cnn};
cnvkit.py coverage {input.bam} {input.antitarget_bed} -o {output.antitarget_cnn};
        """

rule create_reference:
    input:
        cnn = expand(cnv_dir + "{sample}.{prefix}coverage.cnn", sample=samples, prefix=["target","antitarget"]),
        ref = reffasta
    output:
        ref_cnn = cnv_dir + config["analysis"]["case_id"] + "_PON_reference.cnn",
        txt = cnv_dir + "PON." + "reference" + ".done"
    singularity:
        Path(singularity_image, "varcall_cnvkit.sif").as_posix()
    benchmark:
        Path(benchmark_dir, "cnvkit.reference.tsv").as_posix()
    shell:
        """
source activate varcall_cnvkit;
cnvkit.py reference {input.cnn} --fasta {input.ref} -o {output.ref_cnn} && touch {output.txt} ;
        """
