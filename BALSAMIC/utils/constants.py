"""This file contains constants variables used by BALSAMIC"""
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
ANALYSIS_TYPES = ["paired", "single", "umi", "qc", "pon"]
WORKFLOW_SOLUTION = ["BALSAMIC", "Sentieon", "DRAGEN", "Sentieon_umi"]
SEQUENCING_TYPE = ["wgs", "targeted"]

# Variantcaller parameters
VARCALL_PARAMS = {
    "tnscope": {
        "tumor": "--min_init_tumor_lod 1.0 --min_tumor_lod 8",
        "normal": "--min_init_normal_lod 0.5 --min_normal_lod 1.0",
    }
}
# Configuration of VCF settings
VCF_DICT = {
    "TNscope_umi": {
        "mutation": "somatic",
        "type": "SNV",
        "analysis_type": ["single", "paired"],
        "sequencing_type": ["targeted"],
        "workflow_solution": ["Sentieon_umi"],
    },
    "tnscope": {
        "mutation": "somatic",
        "type": "SNV",
        "analysis_type": ["paired", "single"],
        "sequencing_type": ["targeted", "wgs"],
        "workflow_solution": ["Sentieon"],
    },
    "tnhaplotyper": {
        "mutation": "somatic",
        "type": "SNV",
        "analysis_type": ["paired", "single"],
        "sequencing_type": ["targeted", "wgs"],
        "workflow_solution": ["Sentieon"],
    },
    "dnascope": {
        "mutation": "germline",
        "type": "SNV",
        "analysis_type": ["paired", "single"],
        "sequencing_type": ["targeted", "wgs"],
        "workflow_solution": ["Sentieon"],
    },
    "manta": {
        "mutation": "somatic",
        "type": "SV",
        "analysis_type": ["paired", "single"],
        "sequencing_type": ["targeted", "wgs"],
        "workflow_solution": ["BALSAMIC"],
    },
    "cnvkit": {
        "mutation": "somatic",
        "type": "CNV",
        "analysis_type": ["paired", "single"],
        "sequencing_type": ["targeted"],
        "workflow_solution": ["BALSAMIC"],
    },
    "vardict": {
        "mutation": "somatic",
        "type": "SNV",
        "analysis_type": ["paired", "single"],
        "sequencing_type": ["targeted"],
        "workflow_solution": ["BALSAMIC"],
    },
    "manta_germline": {
        "mutation": "germline",
        "type": "SV",
        "analysis_type": ["paired", "single"],
        "sequencing_type": ["targeted", "wgs"],
        "workflow_solution": ["BALSAMIC"],
    },
    "haplotypecaller": {
        "mutation": "germline",
        "type": "SNV",
        "analysis_type": ["paired", "single"],
        "sequencing_type": ["targeted"],
        "workflow_solution": ["BALSAMIC"],
    },
    "delly": {
        "mutation": "somatic",
        "type": "SV",
        "analysis_type": ["paired", "single"],
        "sequencing_type": ["wgs", "targeted"],
        "workflow_solution": ["BALSAMIC"],
    },
    "ascat": {
        "mutation": "somatic",
        "type": "CNV",
        "analysis_type": ["paired"],
        "sequencing_type": ["wgs"],
        "workflow_solution": ["BALSAMIC"],
    },
}

# Minimum required QC-values from HS metrics to be able to pass analysis
HSMETRICS_QC_CHECK = {
    "gicfdna_3.1_hg19_design.bed": {
        "MEAN_TARGET_COVERAGE": 500,
        "FOLD_80_BASE_PENALTY": 1.5,
        "PCT_OFF_BAIT": 0.35,
    },
    "gmcksolid_4.1_hg19_design.bed": {
        "MEAN_TARGET_COVERAGE": 500,
        "FOLD_80_BASE_PENALTY": 1.7,
        "PCT_OFF_BAIT": 0.3,
    },
    "gmsmyeloid_5.2_hg19_design.bed": {
        "MEAN_TARGET_COVERAGE": 1000,
        "FOLD_80_BASE_PENALTY": 1.5,
        "PCT_OFF_BAIT": 0.4,
    },
    "lymphoma_6.1_hg19_design.bed": {
        "MEAN_TARGET_COVERAGE": 1000,
        "FOLD_80_BASE_PENALTY": 1.5,
        "PCT_OFF_BAIT": 0.35,
    },
    "gmslymphoid_7.1_hg19_design.bed": {
        "MEAN_TARGET_COVERAGE": 1000,
        "FOLD_80_BASE_PENALTY": 1.5,
        "PCT_OFF_BAIT": 0.35,
    },
    "twistexomerefseq_9.1_hg19_design.bed": {
        "MEAN_TARGET_COVERAGE": 100,
        "FOLD_80_BASE_PENALTY": 1.8,
        "PCT_OFF_BAIT": 0.25,
    },
    "wgs": {"MEAN_TARGET_COVERAGE": 30},
    "METRIC_CRITERIA": {
        "MEAN_TARGET_COVERAGE": "gt",
        "FOLD_80_BASE_PENALTY": "lt",
        "PCT_OFF_BAIT": "lt",
    },
}

# Configuration of VARDICT settings

VARDICT_SETTINGS = {
    "AD": {"tag_value": 5, "filter_name": "balsamic_low_tumor_ad", "field": "INFO"},
    "DP": {
        "tag_value": 100,
        "filter_name": "balsamic_low_tumor_dp",
        "field": "INFO",
    },
    "MQ": {"tag_value": 40, "filter_name": "balsamic_low_mq", "field": "INFO"},
    "AF_max": {"tag_value": 1, "filter_name": "balsamic_af_one", "field": "INFO"},
    "AF_min": {"tag_value": 0.01, "filter_name": "balsamic_low_af", "field": "INFO"},
    "pop_freq": {
        "tag_value": 0.005,
        "filter_name": "balsamic_high_pop_freq",
        "field": "INFO",
    },
    "varcaller_name": "VarDict",
    "filter_type": "general",
    "analysis_type": "tumor_only",
    "description": "General purpose filters used for filtering VarDict",
}

# Configuration for SENTIEON settings:

SENTIEON_VARCALL_SETTINGS = {
    "AD": {"tag_value": 3, "filter_name": "balsamic_low_tumor_ad", "field": "FORMAT"},
    "DP": {
        "tag_value": 10,
        "filter_name": "balsamic_low_tumor_dp",
        "field": "FORMAT",
    },
    "AF_max": {"tag_value": 1, "filter_name": "balsamic_af_one", "field": "FORMAT"},
    "AF_min": {"tag_value": 0.05, "filter_name": "balsamic_low_af", "field": "FORMAT"},
    "pop_freq": {
        "tag_value": 0.001,
        "filter_name": "balsamic_high_pop_freq",
        "field": "INFO",
    },
    "qss": {
        "tag_value": 20,
        "filter_name": "balsamic_low_quality_scores",
        "field": "FORMAT",
    },
    "strand_reads": {
        "tag_value": 0,
        "filter_name": "balsamic_low_strand_read_counts",
        "field": "FORMAT",
    },
    "sor": {
        "tag_value": 3,
        "filter_name": "balsamic_high_strand_oddsratio",
        "field": "INFO",
    },
    "varcaller_name": "sentieon",
    "filter_type": "general",
    "analysis_type": "tumor_only",
    "description": "General purpose filters used for filtering tnscope and tnhaplotyper",
}

# reference related constants
VALID_REF_FORMAT = ["fasta", "vcf", "text", "gtf", "gff"]
VALID_GENOME_VER = ["hg19", "hg38"]

# reference files
REFERENCE_FILES = {
    "hg38": {
        "reference_genome": {
            "url": "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta",
            "file_type": "fasta",
            "gzip": False,
            "genome_version": "hg38",
            "output_file": "Homo_sapiens_assembly38.fasta",
            "output_path": "genome",
        },
        "dbsnp": {
            "url": "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf",
            "file_type": "vcf",
            "gzip": False,
            "genome_version": "hg38",
            "output_file": "Homo_sapiens_assembly38.dbsnp138.vcf",
            "output_path": "variants",
        },
        "hc_vcf_1kg": {
            "url": "gs://genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
            "file_type": "vcf",
            "gzip": True,
            "genome_version": "hg38",
            "output_file": "1000G_phase1.snps.high_confidence.hg38.vcf",
            "output_path": "variants",
        },
        "mills_1kg": {
            "url": "gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
            "file_type": "vcf",
            "gzip": True,
            "genome_version": "hg38",
            "output_file": "Mills_and_1000G_gold_standard.indels.hg38.vcf",
            "output_path": "variants",
        },
        "known_indel_1kg": {
            "url": "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz",
            "file_type": "vcf",
            "gzip": True,
            "genome_version": "hg38",
            "output_file": "Homo_sapiens_assembly38.known_indels.vcf",
            "output_path": "variants",
        },
        "vcf_1kg": {
            "url": "gs://genomics-public-data/resources/broad/hg38/v0/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf",
            "file_type": "vcf",
            "gzip": False,
            "genome_version": "hg38",
            "output_file": "1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf",
            "output_path": "variants",
        },
        "gnomad_variant": {
            "url": "gs://gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz",
            "file_type": "vcf",
            "gzip": False,
            "genome_version": "hg38",
            "output_file": "gnomad.genomes.r2.1.1.sites.vcf.bgz",
            "output_path": "variants",
        },
        "gnomad_variant_index": {
            "url": "gs://gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz.tbi",
            "file_type": "vcf",
            "gzip": False,
            "genome_version": "hg38",
            "output_file": "gnomad.genomes.r2.1.1.sites.vcf.bgz.tbi",
            "output_path": "variants",
        },
        "cosmicdb": {
            "url": "https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v94/VCF/CosmicCodingMuts.vcf.gz",
            "file_type": "vcf",
            "gzip": True,
            "genome_version": "hg38",
            "output_file": "cosmic_coding_muts_v94.vcf",
            "output_path": "variants",
        },
        "wgs_calling": {
            "url": "gs://genomics-public-data/resources/broad/hg38/v0/wgs_calling_regions.hg38.interval_list",
            "file_type": "text",
            "gzip": False,
            "genome_version": "hg38",
            "output_file": "wgs_calling_regions.v1",
            "output_path": "genome",
        },
        "genome_chrom_size": {
            "url": "https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes",
            "file_type": "text",
            "gzip": False,
            "genome_version": "hg38",
            "output_file": "hg38.chrom.sizes",
            "output_path": "genome",
        },
        "refgene_txt": {
            "url": "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz",
            "file_type": "text",
            "gzip": True,
            "genome_version": "hg38",
            "output_file": "refGene.txt",
            "output_path": "genome",
        },
        "refgene_sql": {
            "url": "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.sql",
            "file_type": "text",
            "gzip": False,
            "genome_version": "hg38",
            "output_file": "refGene.sql",
            "output_path": "genome",
        },
        "rankscore": {
            "url": "https://raw.githubusercontent.com/Clinical-Genomics/reference-files/master/cancer/rank_model/cancer_rank_model_-v0.1-.ini",
            "file_type": "text",
            "gzip": False,
            "genome_version": "hg38",
            "output_file": "cancer_rank_model_-v0.1-.ini",
            "output_path": "genome",
        },
        "access_regions": {
            "url": "https://raw.githubusercontent.com/etal/cnvkit/master/data/access-5k-mappable.hg19.bed",
            "file_type": "text",
            "gzip": False,
            "genome_version": "hg38",
            "output_file": "access_5kb_hg38.txt",
            "output_path": "genome",
        },
        "delly_exclusion": {
            "url": "https://raw.githubusercontent.com/dellytools/delly/master/excludeTemplates/human.hg38.excl.tsv",
            "file_type": "text",
            "gzip": False,
            "genome_version": "hg38",
            "output_file": "delly_exclusion.tsv",
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
            "url": "gs://gatk-legacy-bundles/b37/1000G_phase1.snps.high_confidence.b37.vcf.gz",
            "file_type": "vcf",
            "gzip": True,
            "genome_version": "hg19",
            "output_file": "1kg_phase1_snps_high_confidence_b37.vcf",
            "output_path": "variants",
        },
        "mills_1kg": {
            "url": "gs://gatk-legacy-bundles/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz",
            "file_type": "vcf",
            "gzip": True,
            "genome_version": "hg19",
            "output_file": "mills_1kg_index.vcf",
            "output_path": "variants",
        },
        "known_indel_1kg": {
            "url": "gs://gatk-legacy-bundles/b37/1000G_phase1.indels.b37.vcf.gz",
            "file_type": "vcf",
            "gzip": True,
            "genome_version": "hg19",
            "output_file": "1kg_known_indels_b37.vcf",
            "output_path": "variants",
        },
        "vcf_1kg": {
            "url": "gs://genomics-public-data/ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz",
            "file_type": "vcf",
            "gzip": True,
            "genome_version": "hg19",
            "output_file": "1k_genome_wgs_p1_v3_all_sites.vcf",
            "output_path": "variants",
        },
        "gnomad_variant": {
            "url": "gs://gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.vcf.bgz",
            "file_type": "vcf",
            "gzip": False,
            "genome_version": "hg19",
            "output_file": "gnomad.genomes.r2.1.1.sites.vcf.bgz",
            "output_path": "variants",
        },
        "gnomad_variant_index": {
            "url": "gs://gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.vcf.bgz.tbi",
            "file_type": "vcf",
            "gzip": False,
            "genome_version": "hg19",
            "output_file": "gnomad.genomes.r2.1.1.sites.vcf.bgz.tbi",
            "output_path": "variants",
        },
        "cosmicdb": {
            "url": "https://cancer.sanger.ac.uk/cosmic/file_download/GRCh37/cosmic/v94/VCF/CosmicCodingMuts.vcf.gz",
            "file_type": "vcf",
            "gzip": True,
            "genome_version": "hg19",
            "output_file": "cosmic_coding_muts_v94.vcf",
            "output_path": "variants",
        },
        "wgs_calling": {
            "url": "gs://gatk-legacy-bundles/b37/wgs_calling_regions.v1.interval_list",
            "file_type": "text",
            "gzip": False,
            "genome_version": "hg19",
            "output_file": "wgs_calling_regions.v1",
            "output_path": "genome",
        },
        "genome_chrom_size": {
            "url": "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes",
            "file_type": "text",
            "gzip": False,
            "genome_version": "hg19",
            "output_file": "hg19.chrom.sizes",
            "output_path": "genome",
        },
        "refgene_txt": {
            "url": "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz",
            "file_type": "text",
            "gzip": True,
            "genome_version": "hg19",
            "output_file": "refGene.txt",
            "output_path": "genome",
        },
        "refgene_sql": {
            "url": "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.sql",
            "file_type": "text",
            "gzip": False,
            "genome_version": "hg19",
            "output_file": "refGene.sql",
            "output_path": "genome",
        },
        "rankscore": {
            "url": "https://raw.githubusercontent.com/Clinical-Genomics/reference-files/master/cancer/rank_model/cancer_rank_model_-v0.1-.ini",
            "file_type": "text",
            "gzip": False,
            "genome_version": "hg19",
            "output_file": "cancer_rank_model_-v0.1-.ini",
            "output_path": "genome",
        },
        "access_regions": {
            "url": "https://raw.githubusercontent.com/etal/cnvkit/master/data/access-5k-mappable.hg19.bed",
            "file_type": "text",
            "gzip": False,
            "genome_version": "hg19",
            "output_file": "access_5kb_hg19.txt",
            "output_path": "genome",
        },
        "delly_exclusion": {
            "url": "https://raw.githubusercontent.com/dellytools/delly/master/excludeTemplates/human.hg19.excl.tsv",
            "file_type": "text",
            "gzip": False,
            "genome_version": "hg19",
            "output_file": "delly_exclusion.tsv",
            "output_path": "genome",
        },
        "ascat_gccorrection": {
            "url": "https://raw.githubusercontent.com/Clinical-Genomics/reference-files/12a6c760fd542c02de2cda286b6245e46f4b6a97/cancer/references/GRCh37_SnpGcCorrections.tsv.gz",
            "file_type": "text",
            "gzip": True,
            "genome_version": "hg19",
            "output_file": "GRCh37_SnpGcCorrections.tsv",
            "output_path": "genome",
        },
        "ascat_chryloci": {
            "url": "https://raw.githubusercontent.com/Clinical-Genomics/reference-files/12a6c760fd542c02de2cda286b6245e46f4b6a97/cancer/references/GRCh37_d5_Y.loci",
            "file_type": "text",
            "gzip": False,
            "genome_version": "hg19",
            "output_file": "GRCh37_Y.loci",
            "output_path": "genome",
        },
    },
}

workflow_params = {
    "common": {
        "pcr_model": "NONE",
        "align_header": "'@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINAi'",
        "min_mapq": "20",
        "picard_RG_normal": " ".join(
            [
                "RGPU=ILLUMINAi",
                "RGID=NORMAL",
                "RGSM=NORMAL",
                "RGPL=ILLUMINAi",
                "RGLB=ILLUMINAi",
            ]
        ),
        "picard_RG_tumor": " ".join(
            [
                "RGPU=ILLUMINAi",
                "RGID=TUMOR",
                "RGSM=TUMOR",
                "RGPL=ILLUMINAi",
                "RGLB=ILLUMINAi",
            ]
        ),
    },
    "vardict": {
        "allelic_frequency": "0.001",
        "max_pval": "0.9",
        "max_mm": "4.5",
        "column_info": "-c 1 -S 2 -E 3 -g 4",
    },
    "vep": {
        "vep_filters": "--compress_output bgzip --vcf --everything --allow_non_variant --dont_skip --buffer_size 10000 --format vcf --offline --variant_class --merged --cache --verbose --force_overwrite"
    },
    "umicommon": {
        "align_header": "'@RG\\tID:{sample}\\tSM:{sample}\\tLB:TargetPanel\\tPL:ILLUMINA'",
        "align_intbases": 1000000,
        "filter_tumor_af": 0.0005,
    },
    "umiconsensuscall": {
        "align_format": "BAM",
        "filter_minreads": "3,1,1",
        "tag": "XR",
    },
    "umiextract": {"read_structure": "-d '3M2S+T,3M2S+T'"},
    "tnscope_umi": {
        "algo": "TNscope",
        "min_tumorLOD": 4,
        "init_tumorLOD": 0.5,
        "error_rate": 5,
        "prunefactor": 3,
        "disable_detect": "sv",
    },
}

# list of bioinfo tools for each conda env
VALID_CONTAINER_CONDA_NAME = {
    "align_qc",
    "annotate",
    "coverage_qc",
    "varcall_py36",
    "varcall_py27",
    "varcall_cnvkit",
    "varcall_delly",
    "ascatngs",
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
    "strelka": "varcall_py27",
    "manta": "varcall_py27",
    "cnvkit": "varcall_cnvkit",
    "delly": "varcall_delly",
    "ascatNgs": "ascatngs",
    "sentieon": "sentieon",
}

REPORT_MODEL = {
    "qc": {
        "MEDIAN_TARGET_COVERAGE": {
            "sv": "Mediansekvensdjup [x]",
            "en": "Median sequencing depth [x]",
            "decimal": 0,
        },
        "FOLD_80_BASE_PENALTY": {
            "sv": "Fold 80 base penalty",
            "en": "Fold 80 base penalty",
            "decimal": 2,
        },
        "MEAN_INSERT_SIZE": {
            "sv": "Fragmentlängd, medel [baspar]",
            "en": "Mean insert size [base pair]",
            "decimal": 2,
        },
    },
    "coverage": {
        "PCT_TARGET_BASES_50X": {
            "sv": "Täckningsgrad [50X]",
            "en": "Target coverage [50X]",
            "decimal": 2,
            "as_percent": True,
        },
        "PCT_TARGET_BASES_100X": {
            "sv": "Täckningsgrad [100X]",
            "en": "Target coverage [100X]",
            "decimal": 2,
            "as_percent": True,
        },
        "PCT_TARGET_BASES_250X": {
            "sv": "Täckningsgrad [250X]",
            "en": "Target coverage [250X]",
            "decimal": 2,
            "as_percent": True,
        },
        "PCT_TARGET_BASES_500X": {
            "sv": "Täckningsgrad [500X]",
            "en": "Target coverage [500X]",
            "decimal": 2,
            "as_percent": True,
        },
        "PCT_TARGET_BASES_1000X": {
            "sv": "Täckningsgrad [1000X]",
            "en": "Target coverage [1000X]",
            "decimal": 2,
            "as_percent": True,
        },
    },
}
