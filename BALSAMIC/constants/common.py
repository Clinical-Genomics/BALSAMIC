"""This file contains constants variables used by BALSAMIC"""
import operator
import sys
from pathlib import Path

# DOCKER hub path
BALSAMIC_DOCKER_PATH = "docker://clinicalgenomics/balsamic"

# BALSAMIC base dir
BALSAMIC_BASE_DIR = Path(sys.modules["BALSAMIC"].__file__).parent.resolve()

# BALSAMIC scripts dir
BALSAMIC_SCRIPTS = Path(BALSAMIC_BASE_DIR, "assets/scripts").as_posix()

# Path to containers directory containing YAML files for conda installation for each one
CONTAINERS_CONDA_ENV_PATH = Path(BALSAMIC_BASE_DIR / "containers").as_posix()

# Path to rule files to be accessed by Snakemake
RULE_DIRECTORY = BALSAMIC_BASE_DIR.as_posix()

# Sentieon specific
SENTIEON_DNASCOPE = Path(
    BALSAMIC_BASE_DIR
    / "assets/sentieon_models/SentieonDNAscopeModelBeta0.4a-201808.05.model"
).as_posix()
SENTIEON_TNSCOPE = Path(
    BALSAMIC_BASE_DIR
    / "assets/sentieon_models/SentieonTNscopeModel_GiAB_HighAF_LowFP-201711.05.model"
)

# Analysis related constants
GENDER_OPTIONS = ["female", "male"]
ANALYSIS_TYPES = ["paired", "single", "pon"]
ANALYSIS_WORKFLOW = ["balsamic", "balsamic-qc", "balsamic-umi"]
SEQUENCING_TYPE = ["wgs", "targeted"]
SAMPLE_TYPE = ["tumor", "normal"]
MUTATION_CLASS = ["somatic", "germline"]
MUTATION_TYPE = ["SNV", "SV", "CNV"]
WORKFLOW_SOLUTION = ["BALSAMIC", "Sentieon", "DRAGEN", "Sentieon_umi"]

# Data related constants
FASTQ_SUFFIXES = {
    '1': {'fwd': '_1.fastq.gz', 'rev': '_2.fastq.gz'},
    '2': {'fwd': '_R1_001.fastq.gz', 'rev': '_R2_001.fastq.gz'}
}

# list of bioinfo tools for each conda env
VALID_CONTAINER_CONDA_NAME = {
    "align_qc",
    "annotate",
    "coverage_qc",
    "varcall_py3",
    "varcall_py27",
    "varcall_cnvkit",
    "delly",
    "ascatNgs",
    "balsamic",
    "vcf2cytosure",
    "cnvpytor",
    "somalier",
}

BIOINFO_TOOL_ENV = {
    "bedtools": "align_qc",
    "bwa": "align_qc",
    "fastqc": "align_qc",
    "samtools": "align_qc",
    "picard": "align_qc",
    "compress": "align_qc",
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
