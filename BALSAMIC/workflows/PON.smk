# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from pathlib import Path
import glob
import tempfile
import os


def get_fastqs(wildcards):
    return glob.glob(fastq_dir + "/*_" + wildcards.sample  + "*_R_[1-2].fastq.gz")

def picard_flag(picarddup):
    if picarddup == "mrkdup":
        return "FALSE"
    else:
        return "TRUE"

shell.prefix("set -eo pipefail; ")

localrules: all

case_id = config["analysis"]["case_id"]
pon_ids = config["analysis"]["pon_ids"]
fastq_dir = config["analysis"]["fastq_path"]
analysis_dir = config["analysis"]["analysis_dir"]
reffasta = config["reference"]["reference_genome"]
refflat = config["reference"]["refflat"]
access_5kb_hg19 = config["reference"]["access_regions"]
target_bed = config["panel"]["capture_kit"]
singularity_image = config["singularity"]["image"]

bam_dir = os.path.join(analysis_dir, "bam", "")
cnv_dir = os.path.join(analysis_dir, "cnv", "")
tmp_dir = os.path.join(analysis_dir, "tmp", "" )
Path.mkdir(Path(tmp_dir), exist_ok=True)

picarddup="mrkdup"
picard_extra_normal=" ".join(["RGPU=ILLUMINAi", "RGID=PON","RGSM=PON", "RGPL=ILLUMINAi", "RGLB=ILLUMINAi"])

ALL_COVS = expand(cnv_dir + "{sample}.{cov}coverage.cnn", sample=pon_ids, cov=['target','antitarget'])
ALL_REFS = expand(cnv_dir + "{cov}.bed", cov=['target','antitarget'])
ALL_PON = expand(cnv_dir + config["analysis"]["case_id"] + "_PON_reference.cnn")
PON_DONE = expand(cnv_dir + "PON." + "reference" + ".done")
CNV_DONE = expand(cnv_dir + "CNV." + "tumor" + ".done")

rule all:
    input: ALL_REFS + ALL_COVS +  PON_DONE + CNV_DONE

rule align_bwa_mem:
    input:
        fa = reffasta,
        reads = get_fastqs,
        refidx = expand(reffasta + ".{prefix}", prefix=["amb","ann","bwt","pac","sa"])
    output:
        bamout = bam_dir + "{sample}.sorted.bam"
    params:
        bam_header = "'@RG\\tID:" +  "{sample}" + "\\tSM:" + "{sample}" + "\\tPL:ILLUMINAi'",
        conda = "align_qc",
        tmpdir = tempfile.mkdtemp(prefix=tmp_dir),
    threads: 18
    singularity: 
        Path(singularity_image + "align_qc.sif").as_posix()
    shell:
        """
source activate {params.conda};
mkdir -p {params.tmpdir};
export TMPDIR={params.tmpdir};
bwa mem \
-t {threads} \
-R {params.bam_header}  \
-M \
-v 1 \
{input.fa} {input.reads} \
| samtools sort \
 -T {params.tmpdir} \
--threads {threads} \
--output-fmt BAM \
-o {output.bamout} - ;
samtools index -@ {threads} {output.bamout};
rm -rf {params.tmpdir};
        """

rule picard_markduplicates:
    input:
        bam_dir + "{sample}.sorted.bam"
    output:
        mrkdup = bam_dir + "{sample}.sorted." + picarddup  + ".bam",
        stats =  bam_dir + "{sample}.sorted." + picarddup + ".txt"
    params:
        conda = "align_qc",
        tmpdir = tempfile.mkdtemp(prefix=tmp_dir),
        mem = "16g",
        rm_dup = picard_flag(picarddup)
    threads: 12
    singularity: 
        Path(singularity_image + "align_qc.sif").as_posix() 
    shell:
        """
source activate {params.conda};
mkdir -p {params.tmpdir};
export TMPDIR={params.tmpdir};
picard -Djava.io.tmpdir={params.tmpdir} -Xmx{params.mem} \
MarkDuplicates \
INPUT={input} \
OUTPUT={output.mrkdup} \
VALIDATION_STRINGENCY=SILENT \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
REMOVE_DUPLICATES={params.rm_dup} \
METRICS_FILE='{output.stats}';
samtools index {output.mrkdup};
rm -rf {params.tmpdir};
        """


rule create_target:
    input:
        target_bait = target_bed,
        refFlat = refflat,
        access_bed = access_5kb_hg19
    output:
        target_bed = cnv_dir + "target.bed",
        offtarget_bed = cnv_dir + "antitarget.bed"
    singularity:
        Path(singularity_image + "varcall_cnvkit.sif").as_posix()
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
        Path(singularity_image + "varcall_cnvkit.sif").as_posix()
    shell:
        """
source activate varcall_cnvkit;
cnvkit.py coverage {input.bam} {input.target_bed} -o {output.target_cnn};
cnvkit.py coverage {input.bam} {input.antitarget_bed} -o {output.antitarget_cnn};
        """

rule create_reference:
    input:
        cnn = expand(cnv_dir + "{sample}.{prefix}coverage.cnn", sample=pon_ids, prefix=["target","antitarget"]),
        ref = reffasta
    output:
        ref_cnn = cnv_dir + config["analysis"]["case_id"] + "_PON_reference.cnn",
        txt = cnv_dir + "PON." + "reference" + ".done"
    singularity:
        Path(singularity_image + "varcall_cnvkit.sif").as_posix()
    shell:
        """
source activate varcall_cnvkit;
cnvkit.py reference {input.cnn} --fasta {input.ref} -o {output.ref_cnn} && touch {output.txt} ;
        """

## Test rule: for running cnvkit if tumor sample is provided
rule cnvkit_analysis_single:
    input:
        pon =  cnv_dir + config["analysis"]["case_id"] + "_PON_reference.cnn",
        tumor = bam_dir + "tumor.merged.bam",
    output:
        txt = cnv_dir + "CNV." + "tumor" + ".done"
    params:
        cnv_dir = cnv_dir
    singularity:
         Path(singularity_image + "varcall_cnvkit.sif").as_posix()
    shell:
        """
source activate varcall_cnvkit;
cnvkit.py batch {input.tumor} \
--reference {input.pon} \
--output-dir {params.cnv_dir} \
--diagram --scatter && touch {output.txt};
        """
