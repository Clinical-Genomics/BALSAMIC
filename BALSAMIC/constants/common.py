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

# Path to vcfanno toml files
VCFANNO_TOML = Path(
    BALSAMIC_BASE_DIR / "assets" / "vcfanno" / "vcfanno.toml"
).as_posix()

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
MUTATION_CLASS = ["somatic", "germline"]
MUTATION_TYPE = ["SNV", "SV", "CNV"]
ANALYSIS_TYPES = ["paired", "single", "qc", "pon"]
WORKFLOW_SOLUTION = ["BALSAMIC", "Sentieon", "DRAGEN", "Sentieon_umi"]
SEQUENCING_TYPE = ["wgs", "targeted"]


# list of bioinfo tools for each conda env
VALID_CONTAINER_CONDA_NAME = {
    "align_qc",
    "annotate",
    "coverage_qc",
    "varcall_py36",
    "varcall_py27",
    "varcall_cnvkit",
    "delly",
    "ascatNgs",
    "balsamic",
}

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
    "bcftools": "varcall_py36",
    "tabix": "varcall_py36",
    "gatk": "varcall_py36",
    "vardict": "varcall_py36",
    "svdb": "varcall_py36",
    "strelka": "varcall_py27",
    "manta": "varcall_py27",
    "cnvkit": "varcall_cnvkit",
    "delly": "delly",
    "ascatNgs": "ascatNgs",
    "sentieon": "sentieon",
}

VALID_OPS = {
    "lt": operator.lt,
    "le": operator.le,
    "eq": operator.eq,
    "ne": operator.ne,
    "ge": operator.ge,
    "gt": operator.gt,
}
