METRICS = {
    "targeted": {
        "default": {
            "MEAN_INSERT_SIZE": {"condition": None},
            "PERCENT_DUPLICATION": {"condition": None},
            "MEDIAN_TARGET_COVERAGE": {"condition": None},
            "PCT_TARGET_BASES_20X": {"condition": None},
            "PCT_TARGET_BASES_50X": {"condition": None},
            "PCT_TARGET_BASES_100X": {"condition": None},
            "PCT_TARGET_BASES_250X": {"condition": None},
            "PCT_TARGET_BASES_500X": {"condition": None},
            "PCT_TARGET_BASES_1000X": {"condition": None},
            "MEAN_TARGET_COVERAGE": {"condition": None},
            "FOLD_80_BASE_PENALTY": {"condition": None},
            "PCT_OFF_BAIT": {"condition": None},
            "GC_DROPOUT": {"condition": {"norm": "lt", "threshold": 1.00}},
        },
        "gicfdna": {
            "PCT_TARGET_BASES_1000X": {"condition": {"norm": "gt", "threshold": 0.95}},
        },
        "gmcksolid": {
            "PCT_TARGET_BASES_250X": {"condition": {"norm": "gt", "threshold": 0.95}},
        },
        "gmslymphoid": {
            "PCT_TARGET_BASES_500X": {"condition": {"norm": "gt", "threshold": 0.95}},
        },
        "gmsmyeloid": {
            "PCT_TARGET_BASES_500X": {"condition": {"norm": "gt", "threshold": 0.95}},
        },
        "hdcfdna": {
            "PCT_TARGET_BASES_1000X": {"condition": {"norm": "gt", "threshold": 0.95}},
        },
        "gmssolid": {
            "PCT_TARGET_BASES_250X": {"condition": {"norm": "gt", "threshold": 0.95}},
        },
        "lymphoma": {
            "PCT_TARGET_BASES_500X": {"condition": {"norm": "gt", "threshold": 0.90}},
        },
        "lymphomatic": {
            "PCT_TARGET_BASES_500X": {"condition": {"norm": "gt", "threshold": 0.90}},
        },
        "probio": {
            "PCT_TARGET_BASES_250X": {"condition": {"norm": "gt", "threshold": 0.95}},
        },
        "swep53cfDNA": {
            "PCT_TARGET_BASES_1000X": {"condition": {"norm": "gt", "threshold": 0.95}},
        },
        "twistexome": {
            "PCT_TARGET_BASES_20X": {"condition": {"norm": "gt", "threshold": 0.95}},
        },
        "twistexomerefseq": {
            "PCT_TARGET_BASES_20X": {"condition": {"norm": "gt", "threshold": 0.95}},
        },
        "twistexomecomprehensive": {
            "PCT_TARGET_BASES_20X": {"condition": {"norm": "gt", "threshold": 0.95}},
        },
    },
    "wgs": {
        "MEAN_INSERT_SIZE": {"condition": None},
        "MEDIAN_COVERAGE": {"condition": None},
        "FastQC_mqc-generalstats-fastqc-percent_duplicates": {"condition": None},
        "PCT_15X": {"condition": {"norm": "gt", "threshold": 0.95}},  # Normal sample
        "PCT_30X": {"condition": None},
        "PCT_60X": {"condition": {"norm": "gt", "threshold": 0.80}},  # Tumor sample
        "PCT_100X": {"condition": None},
        "FOLD_80_BASE_PENALTY": {"condition": None},
        "PCT_PF_READS_IMPROPER_PAIRS": {"condition": {"norm": "le", "threshold": 0.1}},
    },
    "variants": {
        "NUMBER_OF_SITES": {"condition": {"norm": "lt", "threshold": 50000}},
    },
    "paired": {
        "RELATEDNESS": {"condition": {"norm": "gt", "threshold": 0.80}},
    },
}
