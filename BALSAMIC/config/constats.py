import BALSAMIC
import sys

from pathlib import Path


QC = {
    "picard_rmdup": False,
    "adapter": "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    "min_seq_length": "25",
    "quality_trim": False,
    "adapter_trim": False,
    "umi_trim": False,
    "umi_trim_length": "5"
}

VCF = {
    "manta": {
        "mutation": "somatic",
        "type": "SV"
    },
    "cnvkit": {
        "mutation": "somatic",
        "type": "CNV"
    },
    "vardict": {
        "mutation": "somatic",
        "type": "SNV"
    },
    "pindel": {
        "mutation": "somatic",
        "type": "SV"
    },
    "strelka": {
        "mutation": "somatic",
        "type": "SNV"
    },
    "mutect": {
        "mutation": "somatic",
        "type": "SNV"
    },
    "tnscope": {
        "mutation": "somatic",
        "type": "SNV"
    },
    "tnsnv": {
        "mutation": "somatic",
        "type": "SNV"
    },
    "tnhaplotyper": {
        "mutation": "somatic",
        "type": "SNV"
    },
    "dnascope": {
        "mutation": "germline",
        "type": "SNV"
    },
    "manta_germline": {
        "mutation": "germline",
        "type": "SV"
    },
    "haplotypecaller": {
        "mutation": "germline",
        "type": "SNV"
    },
    "strelka_germline": {
        "mutation": "germline",
        "type": "SNV"
    },
    "vcfmerge": {
        "mutation": "somatic",
        "type": "SNV"
    }
}


ANALYSIS = {
    "case_id": "CASE_BASE_NAME",
    "analysis_type": "ANALYSIS_TYPE[paired, single]",
    "analysis_dir": "BASE_DIR_RESULTS",
    "fastq_path": "BASE_PATH_TO_FASTQ",
    "script": "scripts/",
    "log": "logs/",
    "result": "analysis/",
    "benchmark": "benchmarks/",
    "BALSAMIC_version" : BALSAMIC.__version__
}


BIOINFO_BASE = [
    "bwa", "bcftools", "cutadapt", "fastqc", "gatk", "manta", "picard",
    "sambamba", "strelka", "samtools", "tabix", "vardic"
]




CONDA_ENV_PATH = Path(sys.modules["BALSAMIC"].__file__).parent.resolve() / "conda"

CONDA_ENV_YAML = Path(sys.modules["BALSAMIC"].__file__).parent.resolve() / "config" / "balsamic_env.yaml"

RULE_DIRECTORY = Path(sys.modules["BALSAMIC"].__file__).parent.resolve())