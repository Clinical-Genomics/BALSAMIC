# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

shell.prefix("set -eo pipefail; ")

from pathlib import Path
import glob

localrules: all

case_ids = ["goodboar","funkite","betterwolf","eagerfeline","chiefweevil"]

raw_folder = "/home/proj/production/cancer/cases"
singularity_image = "/home/proj/stage/cancer/balsamic_cache/7.1.7/containers/"
reffasta = "/home/proj/stage/cancer/balsamic_cache/7.1.7/hg19/genome/human_g1k_v37.fasta"
refflat = "/home/proj/stage/cancer/balsamic_cache/7.1.7/hg19/genome/refGene.flat"
target_bed = "/home/proj/production/cancer/reference/target_capture_bed/production/balsamic/gmcksolid_4.1_hg19_design.bed"
access_5kb_hg19 = "/home/proj/long-term-stage/cancer/UMI_analysis_runs_APJ/cancer_test_datasets/PON/data/access-5k-mappable.hg19.bed"
analysis_dir = "/home/proj/long-term-stage/cancer/UMI_analysis_runs_APJ/cancer_test_datasets/PON/PON_analysis"

bam_dir = analysis_dir + "/bam/"
cnv_dir = analysis_dir + "/cnv/"
tmp_dir= analysis_dir + "/tmp/"

picarddup="mrkdup"
picard_extra_normal=" ".join(["RGPU=ILLUMINAi", "RGID=PON","RGSM=PON", "RGPL=ILLUMINAi", "RGLB=ILLUMINAi"])


def get_fastqs(wildcards):
    return glob.glob(raw_folder+"/"+ wildcards.sample+ "/fastq/"+"*_R_[1-2].fastq.gz")

def get_coverage_fls(wildcards):
    return glob.glob(cnv_dir + "/" + "*coverage.cnn")


def picard_flag(picarddup):
    if picarddup == "mrkdup":
        return "FALSE"
    else:
        return "TRUE"

ALL_BAMS = expand(bam_dir + "{sample}.normal.merged.bam", sample=case_ids)
ALL_COVS = expand(cnv_dir + "{sample}.{cov}coverage.cnn", sample=case_ids, cov=['target','antitarget'])
ALL_REFS = expand(cnv_dir + "{cov}.bed", cov=['target','antitarget'])
ALL_PON = expand(cnv_dir + "PON_reference.cnn")

rule all:
    input: ALL_BAMS + ALL_COVS + ALL_COVS  + ALL_PON

rule align_bwa_mem:
    input:
        fa = reffasta,
        reads = get_fastqs,
        refidx = expand(reffasta + ".{prefix}", prefix=["amb","ann","bwt","pac","sa"])
    output:
        bamout = temp(bam_dir + "{sample}.sorted.bam")
    params:
        bam_header = "'@RG\\tID:" +  "{sample}" + "\\tSM:" + "{sample}" + "\\tPL:ILLUMINAi'",
        conda = "align_qc",
        tmpdir = tmp_dir,
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
| samtools sort -T {params.tmpdir} --threads {threads} --output-fmt BAM -o {output.bamout} - ;
samtools index -@ {threads} {output.bamout};
rm -rf {params.tmpdir};
        """

rule MarkDuplicates:
    input:
        Path(bam_dir, "{sample}.sorted.bam").as_posix()
    output:
        mrkdup = Path(bam_dir, "{sample}.sorted." + picarddup  + ".bam").as_posix(),
        stats = Path(bam_dir, "{sample}.sorted." + picarddup + ".txt").as_posix()
    params:
        conda = "align_qc",
        tmpdir = tmp_dir,
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

rule mergeBam:
    input:
        fasta = reffasta,
        bam = Path(bam_dir, "{sample}.sorted." + picarddup  + ".bam").as_posix()
    output:
        bam = Path(bam_dir + "{sample}.normal.merged.bam").as_posix()
    params:
        conda = "align_qc",
        picard = picard_extra_normal,
    threads: 12
    singularity:
        Path(singularity_image + "align_qc.sif").as_posix() 
    shell:
        """
source activate {params.conda};
picard AddOrReplaceReadGroups {params.picard} INPUT={input.bam} OUTPUT={output.bam};
samtools index {output.bam};
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
        bam = Path(bam_dir + "{sample}.normal.merged.bam").as_posix(),
        target_bed = Path(cnv_dir + "target.bed").as_posix(),
        antitarget_bed = Path(cnv_dir + "antitarget.bed").as_posix()
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
        cnn = get_coverage_fls,
        ref = reffasta
    output:
        ref_cnn = Path(cnv_dir + "PON_reference.cnn").as_posix()
    singularity:
        Path(singularity_image + "varcall_cnvkit.sif").as_posix()
    shell:
        """
source activate varcall_cnvkit;
cnvkit.py reference {input.cnn} --fasta {input.ref} -o {output.ref_cnn};
        """

