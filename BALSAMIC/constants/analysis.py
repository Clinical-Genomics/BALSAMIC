"""Balsamic analysis workflow constants."""
import operator

# Analysis related constants
GENDER_OPTIONS = ["female", "male"]
ANALYSIS_TYPES = ["paired", "single", "pon"]
ANALYSIS_WORKFLOW = ["balsamic", "balsamic-qc", "balsamic-umi"]
SEQUENCING_TYPE = ["wgs", "targeted"]
SAMPLE_TYPE = ["tumor", "normal"]
MUTATION_CLASS = ["somatic", "germline"]
MUTATION_TYPE = ["SNV", "SV", "CNV"]
WORKFLOW_SOLUTION = ["BALSAMIC", "Sentieon", "DRAGEN", "Sentieon_umi"]

BIOINFO_TOOL_ENV = {
    "bedtools": "align_qc",
    "bwa": "align_qc",
    "fastqc": "align_qc",
    "samtools": "align_qc",
    "picard": "align_qc",
    "multiqc": "align_qc",
    "fastp": "align_qc",
    "csvkit": "align_qc",
    "ensembl-vep": "annotate",
    "genmod": "annotate",
    "vcfanno": "annotate",
    "sambamba": "coverage_qc",
    "mosdepth": "coverage_qc",
    "bcftools": "varcall_py3",
    "tabix": "varcall_py3",
    "bgzip": "varcall_py3",
    "gatk": "varcall_py3",
    "vardict": "varcall_py3",
    "svdb": "varcall_py3",
    "tiddit": "varcall_py3",
    "cnvpytor": "cnvpytor",
    "manta": "varcall_py27",
    "cnvkit": "varcall_cnvkit",
    "delly": "delly",
    "ascatNgs": "ascatNgs",
    "sentieon": "sentieon",
    "vcf2cytosure": "vcf2cytosure",
    "somalier": "somalier",
}

VALID_OPS = {
    "lt": operator.lt,
    "le": operator.le,
    "eq": operator.eq,
    "ne": operator.ne,
    "ge": operator.ge,
    "gt": operator.gt,
}
