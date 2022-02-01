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
    "wgs": {"FOLD_80_BASE_PENALTY": {"condition": {"norm": "lt", "threshold": 1.8}}},
}
