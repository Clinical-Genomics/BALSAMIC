COMMON_SETTINGS = {
    "pop_freq": {
        "tag_value": 0.005,
        "filter_name": "balsamic_high_pop_freq",
        "field": "INFO",
    },
    "varcaller_name": "None",
    "filter_type": "general",
    "analysis_type": "tumor_only,tumor_normal",
    "description": "General purpose filters used for filtering any variant caller",
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
    "AF_min": {"tag_value": 0.007, "filter_name": "balsamic_low_af", "field": "INFO"},
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
    "pop_freq_umi": {
        "tag_value": 0.02,
        "filter_name": "balsamic_umi_high_pop_freq",
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
    "varcaller_name": "sentieon",
    "filter_type": "general",
    "analysis_type": "tumor_only",
    "description": "General purpose filters used for filtering tnscope and tnhaplotyper",
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
