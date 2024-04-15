
# Configuration of common SNV filter settings
SNV_BCFTOOLS_SETTINGS_COMMON = {
    "swegen_snv_freq": {
        "tag_value": 0.01,
        "filter_name": "SWEGENAF",
        "field": "INFO",
    },
    "loqusdb_clinical_snv_freq": {
        "tag_value": 0.01,
        "filter_name": "Frq",
        "field": "INFO",
    },
    "high_normal_tumor_af_frac": {
        "tag_value": 0.3,
        "filter_name": "high_normal_tumor_af_frac",
        "field": "FORMAT",
    },
    "sor": {
        "tag_value": 3,
        "filter_name": "balsamic_high_strand_oddsratio",
        "field": "INFO",
    },
    "varcaller_name": "sentieon",
    "filter_type": "general",
    "analysis_type": "tumor_only,tumor_normal",
    "description": "General purpose filters used for filtering TNscope",
}

# Configuration of common TGA SNV filter settings
SNV_BCFTOOLS_SETTINGS_TGA = {
    **SNV_BCFTOOLS_SETTINGS_COMMON,
    "pop_freq": {
        "tag_value": 0.005,
        "filter_name": "balsamic_high_pop_freq",
        "field": "INFO",
    },
    "pop_freq_umi": {
        "tag_value": 0.02,
        "filter_name": "balsamic_umi_high_pop_freq",
        "field": "INFO",
    },
    "AF_min": {"tag_value": 0.007, "filter_name": "balsamic_low_af", "field": "INFO"},
    "AD": {"tag_value": 5, "filter_name": "balsamic_low_tumor_ad", "field": "INFO"},
}

# Configuration of unique TGA SNV filter settings for smaller panels
SNV_BCFTOOLS_SETTINGS_PANEL = {
    **SNV_BCFTOOLS_SETTINGS_TGA,
    "DP": {
        "tag_value": 100,
        "filter_name": "balsamic_low_tumor_dp",
        "field": "INFO",
    },
}

# Configuration of unique TGA SNV filter settings for smaller exomes
SNV_BCFTOOLS_SETTINGS_EXOME = {
    **SNV_BCFTOOLS_SETTINGS_TGA,
    "DP": {
        "tag_value": 20,
        "filter_name": "balsamic_low_tumor_dp",
        "field": "INFO",
    },
}

# Configuration of unique WGS SNV filter settings
SNV_BCFTOOLS_SETTINGS_WGS = {
    **SNV_BCFTOOLS_SETTINGS_COMMON,
    "AD": {"tag_value": 3, "filter_name": "balsamic_low_tumor_ad", "field": "FORMAT"},
    "DP": {
        "tag_value": 10,
        "filter_name": "balsamic_low_tumor_dp",
        "field": "FORMAT",
    },
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
}

# Manta bcftools filters
MANTA_FILTER_SETTINGS = {
    "low_pr_sr_count": {
        "tag_value": 4,
        "filter_name": "low_pr_sr_count",
        "field": "FORMAT",
    },
    "varcaller_name": "Manta",
    "filter_type": "general",
    "analysis_type": "tumor_only,tumor_normal",
    "description": "Bcftools filters to set frequency and minimum read support for SV calls",
}

# Configuration for SVDB settings:
SVDB_FILTER_SETTINGS = {
    "swegen_sv_freq": {
        "tag_value": 0.02,
        "filter_name": "SWEGENAF",
        "field": "INFO",
    },
    "loqusdb_clinical_sv_freq": {
        "tag_value": 0.02,
        "filter_name": "Frq",
        "field": "INFO",
    },
    "varcaller_name": "svdb",
    "filter_type": "general",
    "analysis_type": "tumor_only,tumor_normal",
    "description": "General purpose filters used for filtering svdb merged sv",
}
