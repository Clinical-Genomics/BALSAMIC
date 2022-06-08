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
        "type": "SNV",
        "analysis_type": ["single", "paired"],
        "sequencing_type": ["targeted"],
        "workflow_solution": ["Sentieon_umi"],
    },
    "tnscope": {
        "mutation": "somatic",
        "type": "SNV",
        "analysis_type": ["paired", "single"],
        "sequencing_type": ["wgs"],
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
    "dellysv": {
        "mutation": "somatic",
        "type": "SV",
        "analysis_type": ["paired", "single"],
        "sequencing_type": ["targeted", "wgs"],
        "workflow_solution": ["BALSAMIC"],
    },
    "dellycnv": {
        "mutation": "somatic",
        "type": "CNV",
        "analysis_type": ["single"],
        "sequencing_type": ["targeted", "wgs"],
        "workflow_solution": ["BALSAMIC"],
    },
    "ascat": {
        "mutation": "somatic",
        "type": "CNV",
        "analysis_type": ["paired"],
        "sequencing_type": ["wgs"],
        "workflow_solution": ["BALSAMIC"],
    },
    "svdb": {
        "mutation": "somatic",
        "type": "SV",
        "analysis_type": ["paired", "single"],
        "sequencing_type": ["targeted", "wgs"],
        "workflow_solution": ["BALSAMIC"],
    },
}

WORKFLOW_PARAMS = {
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
        "vep_filters": "--compress_output bgzip --vcf --everything --hgvsg --allow_non_variant --dont_skip --buffer_size 30000 --format vcf --offline --variant_class --merged --cache --verbose --force_overwrite"
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
        "padding": 100,
        "disable_detect": "sv",
    },
}
