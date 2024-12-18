# OPTIONAL SOFT FILTERS
MATCHED_NORMAL_FILTER_NAMES = {
    "matched_normal_filter_names": ["germline_risk", "high_normal_tumor_af_frac"]
}

# CALLER SPECIFIC HARD-FILTERS
# ---------------------------------

TNSCOPE_HARDFILTERS = [
    {
        "filter_name": "t_lod_fstar",
        "Description": "Tumor does not meet likelihood threshold",
    },
    {
        "filter_name": "low_t_alt_frac",
        "Description": "Site filtered due to low alt allele fraction",
    },
]

TNSCOPE_TO_HARDFILTERS = TNSCOPE_HARDFILTERS

TNSCOPE_TN_HARDFILTERS = TNSCOPE_HARDFILTERS + [
    {
        "filter_name": "germline_risk",
        "Description": "Evidence indicates this site is germline, not somatic",
    },
]
VARDICT_HARDFILTERS = [
    {"filter_name": "Cluster0bp", "Description": "Two variants are within 0 bp"},
    {
        "filter_name": "InGap",
        "Description": "The variant is in the deletion gap, thus likely false positive",
    },
    {
        "filter_name": "InIns",
        "Description": "The variant is adjacent to an insertion variant",
    },
    {
        "filter_name": "LongMSI",
        "Description": "The somatic variant is flanked by long A/T (>=14)",
    },
    {
        "filter_name": "MSI12",
        "Description": "Variant in MSI region with 12 non-monomer MSI or 13 monomer MSI",
    },
    {
        "filter_name": "NM4.5",
        "Description": "Mean mismatches in reads >= 4.5, thus likely false positive",
    },
    {"filter_name": "Q10", "Description": "Mean Mapping Quality Below 10"},
    {"filter_name": "SN1.5", "Description": "Signal to Noise Less than 1.5"},
    {"filter_name": "d3", "Description": "Total Depth < 3"},
    {"filter_name": "f0.001", "Description": "Allele frequency < 0.001"},
    {"filter_name": "p8", "Description": "Mean Position in Reads Less than 8"},
    {"filter_name": "pSTD", "Description": "Position in Reads has STD of 0"},
    {"filter_name": "q22.5", "Description": "Mean Base Quality Below 22.5"},
    {"filter_name": "v2", "Description": "Var Depth < 2"},
]

VARDICT_TO_HARDFILTERS = VARDICT_HARDFILTERS + [{"filter_name": "AMPBIAS", "Description": "Indicate the variant has amplicon bias"}]

VARDICT_TN_HARDFILTERS = VARDICT_HARDFILTERS

# RESEARCH AND CLINICAL FILTERS
# ---------------------------------------------------

# Configurations for common clinical filters
SNV_BCFTOOOLS_CLINICAL_COMMON = {
    "loqusdb_clinical_snv_freq": {
        "tag_value": 0.01,
        "filter_name": "Frq",
        "field": "INFO",
    },
    "artefact_snv_freq": {
        "tag_value": 0.1,
        "filter_name": "ArtefactFrq",
        "field": "INFO",
    },
    "varcaller_name": "None",
    "filter_type": "general",
    "analysis_type": "tumor_only,tumor_normal",
    "description": "Clinical database filters used for filtering SNVs",
}

# Configurations for common research filters
SNV_BCFTOOOLS_RESEARCH_COMMON = {
    "swegen_snv_freq": {
        "tag_value": 0.01,
        "filter_name": "SWEGENAF",
        "field": "INFO",
    },
    "varcaller_name": "None",
    "filter_type": "general",
    "analysis_type": "tumor_only,tumor_normal",
    "description": "Research database filters used for filtering SNVs",
}

# Configurations for TGA specific research filters
SNV_BCFTOOOLS_RESEARCH_TGA = {
    **SNV_BCFTOOOLS_RESEARCH_COMMON,
    "pop_freq": {
        "tag_value": 0.005,
        "filter_name": "balsamic_high_pop_freq",
        "field": "INFO",
    },
}

# Configurations for UMI specific research filters
SNV_BCFTOOOLS_RESEARCH_UMI = {
    **SNV_BCFTOOOLS_RESEARCH_COMMON,
    "pop_freq": {
        "tag_value": 0.02,
        "filter_name": "balsamic_high_pop_freq",
        "field": "INFO",
    },
}

# Configurations for WGS specific research filters
SNV_BCFTOOOLS_RESEARCH_WGS = {
    **SNV_BCFTOOOLS_RESEARCH_COMMON,
    "pop_freq": {
        "tag_value": 0.001,
        "filter_name": "balsamic_high_pop_freq",
        "field": "INFO",
    },
}


# SNV QUALITY FILTERS
# ---------------------------------------------------

# Configurations for common tumor normal SNV filters
SNV_BCFTOOLS_QUALITY_TUMOR_NORMAL = {
    "high_normal_tumor_af_frac": {
        "tag_value": 0.3,
        "filter_name": "high_normal_tumor_af_frac",
        "field": "FORMAT",
    },
}

# TGA AND WES
# ---------------------------------------------------

# Configurations for common quality filters for TGA
SNV_BCFTOOLS_QUALITY_COMMON_TGA = {
    "AF_min": {"tag_value": 0.005, "filter_name": "balsamic_low_af", "field": "INFO"},
    "AD": {"tag_value": 5, "filter_name": "balsamic_low_tumor_ad", "field": "INFO"},
}
# Configuration of unique TGA SNV filter settings for smaller panels
SNV_BCFTOOLS_QUALITY_TGA = {
    **SNV_BCFTOOLS_QUALITY_COMMON_TGA,
    "DP": {
        "tag_value": 50,
        "filter_name": "balsamic_low_tumor_dp",
        "field": "INFO",
    },
}

# Configuration of unique TGA SNV filter settings for exomes
SNV_BCFTOOLS_QUALITY_TGA_WES = {
    **SNV_BCFTOOLS_QUALITY_COMMON_TGA,
    "DP": {
        "tag_value": 20,
        "filter_name": "balsamic_low_tumor_dp",
        "field": "INFO",
    },
}


# VARDICT
# --------------------------------------------------
# Configurations for TGA VarDict specific quality filters
SNV_BCFTOOLS_QUALITY_VARDICT_COMMON = {
    "varcaller_name": "vardict",
    "filter_type": "general",
    "analysis_type": "tumor_only,tumor_normal",
    "description": "General purpose filters used for filtering SNVs",
    "MQ": {"tag_value": 30, "filter_name": "balsamic_low_mq", "field": "INFO"},
}

# Configurations for TGA VarDict specific quality filters
SNV_BCFTOOLS_QUALITY_TGA_VARDICT = {
    **SNV_BCFTOOLS_QUALITY_TGA,
    **SNV_BCFTOOLS_QUALITY_VARDICT_COMMON,
}

# Configurations for TGA VarDict tumor only specific quality filters
SNV_BCFTOOLS_QUALITY_TGA_VARDICT_TO = {
    **SNV_BCFTOOLS_QUALITY_TGA_VARDICT,
}

# Configurations for TGA VarDict tumor normal specific quality filters
SNV_BCFTOOLS_QUALITY_TGA_VARDICT_TN = {
    **SNV_BCFTOOLS_QUALITY_TGA_VARDICT,
    **SNV_BCFTOOLS_QUALITY_TUMOR_NORMAL,
}

# Configurations for VarDict WES specific quality filters
SNV_BCFTOOLS_QUALITY_WES_VARDICT = {
    **SNV_BCFTOOLS_QUALITY_TGA_WES,
    **SNV_BCFTOOLS_QUALITY_VARDICT_COMMON,
}

# Configurations for VarDict WES tumor only specific quality filters
SNV_BCFTOOLS_QUALITY_WES_VARDICT_TO = {
    **SNV_BCFTOOLS_QUALITY_WES_VARDICT,
}

# Configurations for VarDict WES tumor normal specific quality filters
SNV_BCFTOOLS_QUALITY_WES_VARDICT_TN = {
    **SNV_BCFTOOLS_QUALITY_WES_VARDICT,
    **SNV_BCFTOOLS_QUALITY_TUMOR_NORMAL,
}

# TNSCOPE
# --------------------------------------------------
# Configurations for TGA TNscope specific quality filters
SNV_BCFTOOLS_QUALITY_TNSCOPE_COMMON = {
    "varcaller_name": "tnscope",
    "filter_type": "general",
    "analysis_type": "tumor_only,tumor_normal",
    "description": "General purpose filters used for filtering SNVs",
}

# Configuration of unique TGA SNV filter settings for small panels
SNV_BCFTOOLS_QUALITY_TGA_TNSCOPE = {
    **SNV_BCFTOOLS_QUALITY_TGA,
    **SNV_BCFTOOLS_QUALITY_TNSCOPE_COMMON,
}

# Configuration of unique TGA SNV filter settings for exomes
SNV_BCFTOOLS_QUALITY_WES_TNSCOPE = {
    **SNV_BCFTOOLS_QUALITY_TGA_WES,
    **SNV_BCFTOOLS_QUALITY_TNSCOPE_COMMON,
}

# Configurations for common TGA TNscope tumor only quality filters
SNV_BCFTOOLS_QUALITY_TGA_TNSCOPE_TO_COMMON = {
    "qss": {
        "tag_value": 20,
        "filter_name": "balsamic_low_quality_scores",
        "field": "FORMAT",
    },
    "sor": {
        "tag_value": 2.7,
        "filter_name": "balsamic_high_strand_oddsratio",
        "field": "INFO",
    },
    "rpa": {
        "tag_value": 8,
        "filter_name": "balsamic_high_repeat_unit_number",
        "field": "INFO",
    },
}

# Configurations for common TGA TNscope tumor normal quality filters
SNV_BCFTOOLS_QUALITY_TGA_TNSCOPE_TN_COMMON = {}

# Configurations for TGA TNscope tumor only quality filters
SNV_BCFTOOLS_QUALITY_TGA_TNSCOPE_TO = {
    **SNV_BCFTOOLS_QUALITY_TGA_TNSCOPE,
    **SNV_BCFTOOLS_QUALITY_TGA_TNSCOPE_TO_COMMON,
}

# Configurations for TGA TNscope tumor normal quality filters
SNV_BCFTOOLS_QUALITY_TGA_TNSCOPE_TN = {
    **SNV_BCFTOOLS_QUALITY_TGA_TNSCOPE,
    **SNV_BCFTOOLS_QUALITY_TUMOR_NORMAL,
    **SNV_BCFTOOLS_QUALITY_TGA_TNSCOPE_TN_COMMON,
}

# Configurations for WES TNscope tumor only quality filters
SNV_BCFTOOLS_QUALITY_WES_TNSCOPE_TO = {
    **SNV_BCFTOOLS_QUALITY_TGA_WES,
    **SNV_BCFTOOLS_QUALITY_WES_TNSCOPE,
    **SNV_BCFTOOLS_QUALITY_TGA_TNSCOPE_TO_COMMON,
}

# Configurations for WES TNscope tumor normal quality filters
SNV_BCFTOOLS_QUALITY_WES_TNSCOPE_TN = {
    **SNV_BCFTOOLS_QUALITY_TGA_WES,
    **SNV_BCFTOOLS_QUALITY_WES_TNSCOPE,
    **SNV_BCFTOOLS_QUALITY_TUMOR_NORMAL,
    **SNV_BCFTOOLS_QUALITY_TGA_TNSCOPE_TN_COMMON,
}

# TNSCOPE UMI
# --------------------------------------------------

# Configurations for UMI TNscope tumor only quality filters
SNV_BCFTOOLS_QUALITY_TGA_TNSCOPE_UMI_TO = {
    **SNV_BCFTOOLS_QUALITY_TNSCOPE_COMMON,
    "variantcaller_filters": TNSCOPE_HARDFILTERS,
}

# Configurations for UMI TNscope tumor normal quality filters
SNV_BCFTOOLS_QUALITY_TGA_TNSCOPE_UMI_TN = {
    **SNV_BCFTOOLS_QUALITY_TNSCOPE_COMMON,
    **SNV_BCFTOOLS_QUALITY_TUMOR_NORMAL,
    "variantcaller_filters": TNSCOPE_TN_HARDFILTERS,
}

# WGS
# ---------------------------------------------------

# Configurations for common quality filters for WGS
SNV_BCFTOOLS_QUALITY_COMMON_WGS = {
    **SNV_BCFTOOLS_QUALITY_TNSCOPE_COMMON,
    "AD": {"tag_value": 3, "filter_name": "balsamic_low_tumor_ad", "field": "FORMAT"},
    "DP": {
        "tag_value": 10,
        "filter_name": "balsamic_low_tumor_dp",
        "field": "FORMAT",
    },
    "AF_min": {"tag_value": 0.05, "filter_name": "balsamic_low_af", "field": "FORMAT"},
}


# Configuration for WGS SNV filters
SNV_BCFTOOLS_QUALITY_WGS_TNSCOPE_COMMON = {}

# Configuration for WGS tumor only SNV filters
SNV_BCFTOOLS_QUALITY_WGS_TNSCOPE_TO = {
    **SNV_BCFTOOLS_QUALITY_COMMON_WGS,
    **SNV_BCFTOOLS_QUALITY_WGS_TNSCOPE_COMMON,
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
}

# Configuration for WGS tumor normal SNV filters
SNV_BCFTOOLS_QUALITY_WGS_TNSCOPE_TN = {
    **SNV_BCFTOOLS_QUALITY_COMMON_WGS,
    **SNV_BCFTOOLS_QUALITY_WGS_TNSCOPE_COMMON,
    **SNV_BCFTOOLS_QUALITY_TUMOR_NORMAL,
}

# All filters for WGS tumor only
SNV_FILTERS_WGS_TO = {
    "tnscope": {
        **SNV_BCFTOOLS_QUALITY_WGS_TNSCOPE_TO,
        "variantcaller_filters": TNSCOPE_TO_HARDFILTERS,
    },
    "research": {
        **SNV_BCFTOOOLS_RESEARCH_WGS,
    },
    "clinical": {
        **SNV_BCFTOOOLS_RESEARCH_WGS,
        **SNV_BCFTOOOLS_CLINICAL_COMMON,
    },
}

# All filters for WGS tumor normal
SNV_FILTERS_WGS_TN = {
    "tnscope": {
        **SNV_BCFTOOLS_QUALITY_WGS_TNSCOPE_TN,
        "variantcaller_filters": TNSCOPE_TN_HARDFILTERS,
    },
    "research": {
        **SNV_BCFTOOOLS_RESEARCH_WGS,
    },
    "clinical": {
        **SNV_BCFTOOOLS_RESEARCH_WGS,
        **SNV_BCFTOOOLS_CLINICAL_COMMON,
    },
}

# All filters for WES tumor only
SNV_FILTERS_TGA_WES_TO = {
    "tnscope": {
        **SNV_BCFTOOLS_QUALITY_WES_TNSCOPE_TO,
        "variantcaller_filters": TNSCOPE_TO_HARDFILTERS,
    },
    "vardict": {
        **SNV_BCFTOOLS_QUALITY_WES_VARDICT_TO,
        "variantcaller_filters": VARDICT_TO_HARDFILTERS,
    },
    "research": {
        **SNV_BCFTOOOLS_RESEARCH_TGA,
    },
    "clinical": {
        **SNV_BCFTOOOLS_RESEARCH_TGA,
        **SNV_BCFTOOOLS_CLINICAL_COMMON,
    },
}

# All filters for WES tumor normal
SNV_FILTERS_TGA_WES_TN = {
    "tnscope": {
        **SNV_BCFTOOLS_QUALITY_WES_TNSCOPE_TN,
        "variantcaller_filters": TNSCOPE_TN_HARDFILTERS,
    },
    "vardict": {
        **SNV_BCFTOOLS_QUALITY_WES_VARDICT_TN,
        "variantcaller_filters": VARDICT_TN_HARDFILTERS,
    },
    "research": {
        **SNV_BCFTOOOLS_RESEARCH_TGA,
    },
    "clinical": {
        **SNV_BCFTOOOLS_RESEARCH_TGA,
        **SNV_BCFTOOOLS_CLINICAL_COMMON,
    },
}

# All filters for TGA tumor only
SNV_FILTERS_TGA_TO = {
    "tnscope": {
        **SNV_BCFTOOLS_QUALITY_TGA_TNSCOPE_TO,
        "variantcaller_filters": TNSCOPE_TO_HARDFILTERS,
    },
    "vardict": {
        **SNV_BCFTOOLS_QUALITY_TGA_VARDICT_TO,
        "variantcaller_filters": VARDICT_TO_HARDFILTERS,
    },
    "research": {
        **SNV_BCFTOOOLS_RESEARCH_TGA,
    },
    "clinical": {
        **SNV_BCFTOOOLS_RESEARCH_TGA,
        **SNV_BCFTOOOLS_CLINICAL_COMMON,
    },
}

# All filters for TGA tumor normal
SNV_FILTERS_TGA_TN = {
    "tnscope": {
        **SNV_BCFTOOLS_QUALITY_TGA_TNSCOPE_TN,
        "variantcaller_filters": TNSCOPE_TN_HARDFILTERS,
    },
    "vardict": {
        **SNV_BCFTOOLS_QUALITY_TGA_VARDICT_TN,
        "variantcaller_filters": VARDICT_TN_HARDFILTERS,
    },
    "research": {
        **SNV_BCFTOOOLS_RESEARCH_TGA,
    },
    "clinical": {
        **SNV_BCFTOOOLS_RESEARCH_TGA,
        **SNV_BCFTOOOLS_CLINICAL_COMMON,
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
