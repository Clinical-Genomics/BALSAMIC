"""This file contains constants variables used by BALSAMIC"""
import sys
from pathlib import Path

import BALSAMIC

# Path to conda folder containing YAML files with verions of software usen un BALSAMIC workflow
CONDA_ENV_PATH = Path(
    Path(sys.modules["BALSAMIC"].__file__).parent.resolve() /
    "conda").as_posix()

# Path to config YAML file to be accessed by Snakemake
CONDA_ENV_YAML = Path(
    Path(sys.modules["BALSAMIC"].__file__).parent.resolve() / "config" /
    "balsamic_env.yaml").as_posix()

# Path to rule files to be accessed by Snakemake
RULE_DIRECTORY = (
    Path(sys.modules["BALSAMIC"].__file__).parent.resolve().as_posix() + "/")

# BALSAMIC version
BALSAMIC_VERSION = BALSAMIC.__version__

# Analysis related constants
MUTATION_CLASS = ["somatic", "germline"]
MUTATION_TYPE = ["SNV", "SV", "CNV"]
ANALYSIS_TYPES = ["paired", "single", "umi", "qc"]
WORKFLOW_SOLUTION = ["BALSAMIC", "Sentieon", "DRAGEN"]

# Configuration of VCF settings
VCF_DICT = {
    "tnsnv": {
        "mutation": "somatic",
        "type": "SNV",
        "analysis_type": ["paired", "single"],
        "workflow_solution": ["Sentieon"]
    },
    "tnscope": {
        "mutation": "somatic",
        "type": "SNV",
        "analysis_type": ["paired", "single"],
        "workflow_solution": ["Sentieon"]
    },
    "tnhaplotyper": {
        "mutation": "somatic",
        "type": "SNV",
        "analysis_type": ["paired", "single"],
        "workflow_solution": ["Sentieon"]
    },
    "dnascope": {
        "mutation": "germline",
        "type": "SNV",
        "analysis_type": ["paired", "single"],
        "workflow_solution": ["Sentieon"]
    },
    "manta": {
        "mutation": "somatic",
        "type": "SV",
        "analysis_type": ["paired", "single"],
        "workflow_solution": ["BALSAMIC"]
    },
    "cnvkit": {
        "mutation": "somatic",
        "type": "CNV",
        "analysis_type": ["paired", "single"],
        "workflow_solution": ["BALSAMIC"]
    },
    "mutect": {
        "mutation": "somatic",
        "type": "SNV",
        "analysis_type": ["paired", "single"],
        "workflow_solution": ["BALSAMIC"]
    },
    "vardict": {
        "mutation": "somatic",
        "type": "SNV",
        "analysis_type": ["paired", "single"],
        "workflow_solution": ["BALSAMIC"]
    },
    "strelka": {
        "mutation": "somatic",
        "type": "SNV",
        "analysis_type": ["paired"],
        "workflow_solution": ["BALSAMIC"]
    },
    "manta_germline": {
        "mutation": "germline",
        "type": "SV",
        "analysis_type": ["paired", "single"],
        "workflow_solution": ["BALSAMIC"]
    },
    "haplotypecaller": {
        "mutation": "germline",
        "type": "SNV",
        "analysis_type": ["paired", "single"],
        "workflow_solution": ["BALSAMIC"]
    },
    "strelka_germline": {
        "mutation": "germline",
        "type": "SNV",
        "analysis_type": ["paired", "single"],
        "workflow_solution": ["BALSAMIC"]
    }
}

# Minimum required QC-values from HS metrics to be able to pass analysis
HSMETRICS_QC_CHECK = {
    "gicfdna_3.1_hg19_design.bed": {
        "MEAN_TARGET_COVERAGE": 500,
        "FOLD_80_BASE_PENALTY": 1.5,
        "PCT_OFF_BAIT": 0.35
    },
    "gmcksolid_4.1_hg19_design.bed": {
        "MEAN_TARGET_COVERAGE": 500,
        "FOLD_80_BASE_PENALTY": 1.7,
        "PCT_OFF_BAIT": 0.3
    },
    "gmsmyeloid_5.2_hg19_design.bed": {
        "MEAN_TARGET_COVERAGE": 1000,
        "FOLD_80_BASE_PENALTY": 1.5,
        "PCT_OFF_BAIT": 0.4
    },
    "lymphoma_6.1_hg19_design.bed": {
        "MEAN_TARGET_COVERAGE": 1000,
        "FOLD_80_BASE_PENALTY": 1.5,
        "PCT_OFF_BAIT": 0.35
    },
    "gmslymphoid_7.1_hg19_design.bed": {
        "MEAN_TARGET_COVERAGE": 1000,
        "FOLD_80_BASE_PENALTY": 1.5,
        "PCT_OFF_BAIT": 0.35
    },
    "twistexomerefseq_9.1_hg19_design.bed": {
        "MEAN_TARGET_COVERAGE": 100,
        "FOLD_80_BASE_PENALTY": 1.8,
        "PCT_OFF_BAIT": 0.25
    },
    "wgs": {
        "MEAN_TARGET_COVERAGE": 30
    },
    "METRIC_CRITERIA": {
        "MEAN_TARGET_COVERAGE": "gt",
        "FOLD_80_BASE_PENALTY": "lt",
        "PCT_OFF_BAIT": "lt"
    }
}

# Configuration of VARDICT settings

VARDICT_SETTINGS = {
    "AD": {
        "tag_value": 5,
        "filter_name": "balsamic_low_tumor_ad",
        "field": "INFO"
    },
    "DP": {
        "tag_value": 100,
        "filter_name": "balsamic_low_tumor_dp",
        "field": "INFO",
    },
    "MQ": {
        "tag_value": 50,
        "filter_name": "balsamic_low_mq",
        "field": "INFO"
    },
    "AF_max": {
        "tag_value": 1,
        "filter_name": "balsamic_af_one",
        "field": "INFO"
    },
    "AF_min": {
        "tag_value": 0.02,
        "filter_name": "balsamic_low_af",
        "field": "INFO"
    },
    "varcaller_name": "VarDict",
    "filter_type": "general",
    "analysis_type": "tumor_only",
    "description": "General purpose filters used for filtering VarDict",
}

# reference related constants
VALID_REF_FORMAT = ["fasta", "vcf", "text", "gtf", "gff"]
VALID_GENOME_VER = ["hg19", "hg38"]

# reference files
REFERENCE_FILES = {
    "hg38": {
        "reference_genome": {
            "url":
            "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta",
            "file_type": "fasta",
            "gzip": False,
            "genome_version": "hg38",
            "output_file": "Homo_sapiens_assembly38.fasta",
            "output_path": "genome",
        },
        "dbsnp": {
            "url":
            "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf",
            "file_type": "vcf",
            "gzip": False,
            "genome_version": "hg38",
            "output_file": "Homo_sapiens_assembly38.dbsnp138.vcf",
            "output_path": "variants",
        },
        "hc_vcf_1kg": {
            "url":
            "gs://genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
            "file_type": "vcf",
            "gzip": True,
            "genome_version": "hg38",
            "output_file": "1000G_phase1.snps.high_confidence.hg38.vcf",
            "output_path": "variants",
        },
        "mills_1kg": {
            "url":
            "gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
            "file_type": "vcf",
            "gzip": True,
            "genome_version": "hg38",
            "output_file": "Mills_and_1000G_gold_standard.indels.hg38.vcf",
            "output_path": "variants",
        },
        "known_indel_1kg": {
            "url":
            "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz",
            "file_type": "vcf",
            "gzip": True,
            "genome_version": "hg38",
            "output_file": "Homo_sapiens_assembly38.known_indels.vcf",
            "output_path": "variants",
        },
        "vcf_1kg": {
            "url":
            "gs://genomics-public-data/resources/broad/hg38/v0/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf",
            "file_type": "vcf",
            "gzip": False,
            "genome_version": "hg38",
            "output_file":
            "1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf",
            "output_path": "variants",
        },
        "cosmicdb": {
            "url":
            "https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v92/VCF/CosmicCodingMuts.vcf.gz",
            "file_type": "vcf",
            "gzip": True,
            "genome_version": "hg38",
            "output_file": "cosmic_coding_muts_v92.vcf",
            "output_path": "variants",
        },
        "wgs_calling": {
            "url":
            "gs://genomics-public-data/resources/broad/hg38/v0/wgs_calling_regions.hg38.interval_list",
            "file_type": "text",
            "gzip": False,
            "genome_version": "hg38",
            "output_file": "wgs_calling_regions.v1",
            "output_path": "genome",
        },
        "genome_chrom_size": {
            "url":
            "https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes",
            "file_type": "text",
            "gzip": False,
            "genome_version": "hg38",
            "output_file": "hg38.chrom.sizes",
            "output_path": "genome",
        },
        "refgene_txt": {
            "url":
            "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz",
            "file_type": "text",
            "gzip": True,
            "genome_version": "hg38",
            "output_file": "refGene.txt",
            "output_path": "genome",
        },
        "refgene_sql": {
            "url":
            "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.sql",
            "file_type": "text",
            "gzip": False,
            "genome_version": "hg38",
            "output_file": "refGene.sql",
            "output_path": "genome",
        },
    },
    "hg19": {
        "reference_genome": {
            "url": "gs://gatk-legacy-bundles/b37/human_g1k_v37.fasta.gz",
            "file_type": "fasta",
            "gzip": True,
            "genome_version": "hg19",
            "output_file": "human_g1k_v37.fasta",
            "output_path": "genome",
        },
        "dbsnp": {
            "url": "gs://gatk-legacy-bundles/b37/dbsnp_138.b37.vcf.gz",
            "file_type": "vcf",
            "gzip": True,
            "genome_version": "hg19",
            "output_file": "dbsnp_grch37_b138.vcf",
            "output_path": "variants",
        },
        "hc_vcf_1kg": {
            "url":
            "gs://gatk-legacy-bundles/b37/1000G_phase1.snps.high_confidence.b37.vcf.gz",
            "file_type": "vcf",
            "gzip": True,
            "genome_version": "hg19",
            "output_file": "1kg_phase1_snps_high_confidence_b37.vcf",
            "output_path": "variants",
        },
        "mills_1kg": {
            "url":
            "gs://gatk-legacy-bundles/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz",
            "file_type": "vcf",
            "gzip": True,
            "genome_version": "hg19",
            "output_file": "mills_1kg_index.vcf",
            "output_path": "variants",
        },
        "known_indel_1kg": {
            "url":
            "gs://gatk-legacy-bundles/b37/1000G_phase1.indels.b37.vcf.gz",
            "file_type": "vcf",
            "gzip": True,
            "genome_version": "hg19",
            "output_file": "1kg_known_indels_b37.vcf.gz",
            "output_path": "variants",
        },
        "vcf_1kg": {
            "url":
            "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz",
            "file_type": "vcf",
            "gzip": True,
            "genome_version": "hg19",
            "output_file": "1k_genome_wgs_p1_v3_all_sites.vcf",
            "output_path": "variants",
        },
        "cosmicdb": {
            "url":
            "https://cancer.sanger.ac.uk/cosmic/file_download/GRCh37/cosmic/v90/VCF/CosmicCodingMuts.vcf.gz",
            "file_type": "vcf",
            "gzip": True,
            "genome_version": "hg19",
            "output_file": "cosmic_coding_muts_v90.vcf",
            "output_path": "variants",
        },
        "wgs_calling": {
            "url":
            "gs://gatk-legacy-bundles/b37/wgs_calling_regions.v1.interval_list",
            "file_type": "text",
            "gzip": False,
            "genome_version": "hg19",
            "output_file": "wgs_calling_regions.v1",
            "output_path": "genome",
        },
        "genome_chrom_size": {
            "url":
            "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes",
            "file_type": "text",
            "gzip": False,
            "genome_version": "hg19",
            "output_file": "hg19.chrom.sizes",
            "output_path": "genome",
        },
        "refgene_txt": {
            "url":
            "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz",
            "file_type": "text",
            "gzip": True,
            "genome_version": "hg19",
            "output_file": "refGene.txt",
            "output_path": "genome",
        },
        "refgene_sql": {
            "url":
            "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.sql",
            "file_type": "text",
            "gzip": False,
            "genome_version": "hg19",
            "output_file": "refGene.sql",
            "output_path": "genome",
        },
    },
}
