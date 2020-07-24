"""This file contains consatant variables used by BALSAMIC"""
import sys
from pathlib import Path

import BALSAMIC

#Path to conda folder containing YAML files with verions of software usen un BALSAMIC workflow
CONDA_ENV_PATH = Path(
    Path(sys.modules["BALSAMIC"].__file__).parent.resolve() /
    "conda").as_posix()

#Path to config YAML file to be accessed by Snakemake
CONDA_ENV_YAML = Path(
    Path(sys.modules["BALSAMIC"].__file__).parent.resolve() / "config" /
    "balsamic_env.yaml").as_posix()

#Path to rule files to be accessed by Snakemake
RULE_DIRECTORY = Path(
    sys.modules["BALSAMIC"].__file__).parent.resolve().as_posix() + "/"

#BALSMIC version
BALSAMIC_VERSION = BALSAMIC.__version__

#Configuration of VCF settings
VCF_DICT = {
    "tnsnv": {
        "mutation": "somatic",
        "type": "SNV"
    },
    "manta": {
        "mutation": "somatic",
        "type": "SV"
    },
    "pindel": {
        "mutation": "somatic",
        "type": "SV"
    },
    "cnvkit": {
        "mutation": "somatic",
        "type": "CNV"
    },
    "mutect": {
        "mutation": "somatic",
        "type": "SNV"
    },
    "vardict": {
        "mutation": "somatic",
        "type": "SNV"
    },
    "strelka": {
        "mutation": "somatic",
        "type": "SNV"
    },
    "tnscope": {
        "mutation": "somatic",
        "type": "SNV"
    },
    "vcfmerge": {
        "mutation": "somatic",
        "type": "SNV"
    },
    "dnascope": {
        "mutation": "germline",
        "type": "SNV"
    },
    "tnhaplotyper": {
        "mutation": "somatic",
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
}

#Configuration of VARDICT settings
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
    "description": "General purpose filters used for filtering VarDict"
}

#reference related constants
VALID_REF_FORMAT = ["fasta", "vcf", "text", "gtf", "gff"]
VALID_GENOME_VER = ["hg19", "hg38"]

#reference files
REFERENCE_FILES = {
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
            "file_type": "fasta",
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
            "output_file": "cosmic_coding_muts_v89.vcf",
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
    }
}
