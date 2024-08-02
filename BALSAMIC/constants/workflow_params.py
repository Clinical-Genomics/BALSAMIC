# Variantcaller parameters
VARCALL_PARAMS = {
    "tnscope": {
        "tumor": "--min_init_tumor_lod 1.0 --min_tumor_lod 8",
        "normal": "--min_init_normal_lod 0.5 --min_normal_lod 1.0",
    }
}

# Configuration of VCF settings
VCF_DICT = {
    "tnscope_umi": {
        "mutation": "somatic",
        "mutation_type": "SNV",
        "analysis_type": ["single", "paired"],
        "sequencing_type": ["targeted"],
        "workflow_solution": ["Sentieon_umi"],
    },
    "tnscope": {
        "mutation": "somatic",
        "mutation_type": "SNV",
        "analysis_type": ["paired", "single"],
        "sequencing_type": ["targeted", "wgs"],
        "workflow_solution": ["Sentieon"],
    },
    "dnascope": {
        "mutation": "germline",
        "mutation_type": "SNV",
        "analysis_type": ["paired", "single"],
        "sequencing_type": ["targeted", "wgs"],
        "workflow_solution": ["Sentieon"],
    },
    "manta": {
        "mutation": "somatic",
        "mutation_type": "SV",
        "analysis_type": ["paired", "single"],
        "sequencing_type": ["targeted", "wgs"],
        "workflow_solution": ["BALSAMIC"],
    },
    "cnvkit": {
        "mutation": "somatic",
        "mutation_type": "CNV",
        "analysis_type": ["paired", "single"],
        "sequencing_type": ["targeted"],
        "workflow_solution": ["BALSAMIC"],
    },
    "manta_germline": {
        "mutation": "germline",
        "mutation_type": "SV",
        "analysis_type": ["paired", "single"],
        "sequencing_type": ["targeted", "wgs"],
        "workflow_solution": ["BALSAMIC"],
    },
    "dellysv": {
        "mutation": "somatic",
        "mutation_type": "SV",
        "analysis_type": ["paired", "single"],
        "sequencing_type": ["targeted", "wgs"],
        "workflow_solution": ["BALSAMIC"],
    },
    "tiddit": {
        "mutation": "somatic",
        "mutation_type": "SV",
        "analysis_type": ["single", "paired"],
        "sequencing_type": ["wgs"],
        "workflow_solution": ["BALSAMIC"],
    },
    "dellycnv": {
        "mutation": "somatic",
        "mutation_type": "CNV",
        "analysis_type": ["single", "paired"],
        "sequencing_type": ["targeted", "wgs"],
        "workflow_solution": ["BALSAMIC"],
    },
    "ascat": {
        "mutation": "somatic",
        "mutation_type": "CNV",
        "analysis_type": ["paired"],
        "sequencing_type": ["wgs"],
        "workflow_solution": ["BALSAMIC"],
    },
    "cnvpytor": {
        "mutation": "somatic",
        "mutation_type": "CNV",
        "analysis_type": ["single"],
        "sequencing_type": ["wgs"],
        "workflow_solution": ["BALSAMIC"],
    },
    "igh_dux4": {
        "mutation": "somatic",
        "mutation_type": "SV",
        "analysis_type": ["single", "paired"],
        "sequencing_type": ["wgs"],
        "workflow_solution": ["BALSAMIC"],
    },
    "svdb": {
        "mutation": "somatic",
        "mutation_type": "SV",
        "analysis_type": ["paired", "single"],
        "sequencing_type": ["targeted", "wgs"],
        "workflow_solution": ["BALSAMIC"],
    },
}

SLEEP_BEFORE_START = 600

WORKFLOW_PARAMS = {
    "bam_post_processing": {
        "manta_max_base_quality": 70,
    },
    "bed_pre_processing": {
        "minimum_region_size": 100,
    },
    "common": {
        "header_per_lane": "'@RG\\tID:{fastq_pattern}\\tSM:{sample_type}\\tPL:ILLUMINAi'",
        "header_per_sample": "'@RG\\tID:{sample}\\tSM:{sample_type}\\tPL:ILLUMINAi'",
        "pcr_model": "NONE",
        "min_mapq": "20",
        "picard_fixmate": " ".join(
            [
                "-ADD_MATE_CIGAR true",
                "-MAX_RECORDS_IN_RAM 10000000",
                "-CREATE_INDEX true",
                "-CREATE_MD5_FILE true",
            ]
        ),
        "picard_RG_normal": " ".join(
            [
                "-RGPU ILLUMINAi",
                "-RGID NORMAL",
                "-RGSM NORMAL",
                "-RGPL ILLUMINAi",
                "-RGLB ILLUMINAi",
                "-MAX_RECORDS_IN_RAM 1000000",
                "-CREATE_INDEX true",
                "-CREATE_MD5_FILE true",
            ]
        ),
        "picard_RG_tumor": " ".join(
            [
                "-RGPU ILLUMINAi",
                "-RGID TUMOR",
                "-RGSM TUMOR",
                "-RGPL ILLUMINAi",
                "-RGLB ILLUMINAi",
                "-MAX_RECORDS_IN_RAM 1000000",
                "-CREATE_INDEX true",
                "-CREATE_MD5_FILE true",
            ]
        ),
    },
    "insert_size_metrics": {
        "min_read_ratio": 0.01,
    },
    "manta": {
        "wgs_settings": "",
        "tga_settings": "--exome",
    },
    "mosdepth": {
        "mapq": 20,
        "samflag": 1796,
        "quantize": "0:1:50:150:",
    },
    "sentieon_wgs_metrics": {
        "min_base_qual": 10,
        "cov_threshold": [50, 100, 150, 200, 250],
    },
    "vep": {
        "vep_filters": "--compress_output bgzip --vcf --everything --hgvsg --allow_non_variant --dont_skip --buffer_size 30000 --max_sv_size 249250621 --format vcf --offline --variant_class --merged --cache --verbose --force_overwrite"
    },
    "umicommon": {"align_intbases": 1000000},
    "umiconsensuscall": {
        "align_format": "BAM",
        "filter_minreads": "3,1,1",
        "tag": "XR",
    },
    "umiextract": {"read_structure": "-d '3M2S+T,3M2S+T'"},
    "tnscope_umi": {
        "algo": "TNscope",
        "filter_tumor_af": 0.0005,
        "pcr_model": "NONE",
        "min_tumorLOD": 4,
        "init_tumorLOD": 0.5,
        "error_rate": 5,
        "prunefactor": 3,
        "padding": 100,
        "disable_detect": "sv",
    },
    "tnscope_tga": {
        "algo": "TNscope",
        "filter_tumor_af": 0.0005,
        "pcr_model": "NONE",
        "min_tumorLOD": 4,
        "init_tumorLOD": 0.5,
        "error_rate": 5,
        "prunefactor": 3,
        "padding": 100,
    },
}
