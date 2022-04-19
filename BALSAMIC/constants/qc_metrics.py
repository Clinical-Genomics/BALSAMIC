METRICS = {
    "targeted": {
        "default": {
            "MEAN_INSERT_SIZE": {"condition": None},
            "PERCENT_DUPLICATION": {"condition": None},
            "MEDIAN_TARGET_COVERAGE": {"condition": {"norm": "gt", "threshold": 500}},
            "PCT_TARGET_BASES_50X": {"condition": None},
            "PCT_TARGET_BASES_100X": {"condition": None},
            "PCT_TARGET_BASES_250X": {"condition": None},
            "PCT_TARGET_BASES_500X": {"condition": None},
            "PCT_TARGET_BASES_1000X": {"condition": None},
            "MEAN_TARGET_COVERAGE": {"condition": None},
            "FOLD_80_BASE_PENALTY": {"condition": {"norm": "lt", "threshold": 1.8}},
            "PCT_OFF_BAIT": {"condition": None},
        },
        "gicfdna": {
            "MEDIAN_TARGET_COVERAGE": {"condition": {"norm": "gt", "threshold": 1000}},
            "FOLD_80_BASE_PENALTY": {"condition": {"norm": "lt", "threshold": 1.6}},
        },
        "gmcksolid": {
            "MEDIAN_TARGET_COVERAGE": {"condition": {"norm": "gt", "threshold": 500}},
            "FOLD_80_BASE_PENALTY": {"condition": {"norm": "lt", "threshold": 1.8}},
        },
        "gmsmyeloid": {
            "MEDIAN_TARGET_COVERAGE": {"condition": {"norm": "gt", "threshold": 1000}},
            "FOLD_80_BASE_PENALTY": {"condition": {"norm": "lt", "threshold": 1.6}},
        },
        "lymphoma": {
            "MEDIAN_TARGET_COVERAGE": {"condition": {"norm": "gt", "threshold": 1000}},
            "FOLD_80_BASE_PENALTY": {"condition": {"norm": "lt", "threshold": 1.6}},
        },
        "gmslymphoid": {
            "MEDIAN_TARGET_COVERAGE": {"condition": {"norm": "gt", "threshold": 1000}},
            "FOLD_80_BASE_PENALTY": {"condition": {"norm": "lt", "threshold": 1.6}},
        },
        "twistexome": {
            "MEDIAN_TARGET_COVERAGE": {"condition": {"norm": "gt", "threshold": 100}},
            "FOLD_80_BASE_PENALTY": {"condition": {"norm": "lt", "threshold": 1.8}},
        },
    },
    "wgs": {
        "MEAN_INSERT_SIZE": {"condition": None},
        "MEDIAN_COVERAGE": {"condition": None},
        "FastQC_mqc-generalstats-fastqc-percent_duplicates": {"condition": None},
        "PCT_15X": {"condition": None},
        "PCT_30X": {"condition": None},
        "PCT_60X": {"condition": None},
        "PCT_100X": {"condition": None},
        "FOLD_80_BASE_PENALTY": {"condition": {"norm": "lt", "threshold": 1.8}},
    },
}
