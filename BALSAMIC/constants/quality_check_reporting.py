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

METRICS = {
    "qc": {
        "targeted": {
            "multiqc_picard_insertSize.json": {
                "MEAN_INSERT_SIZE": {"condition": None},
            },
            "multiqc_picard_dups.json": {"PERCENT_DUPLICATION": {"condition": None}},
            "multiqc_picard_HsMetrics.json": {
                "MEAN_TARGET_COVERAGE": {
                    "condition": {"norm": "gt", "threshold": 500.0}
                },
                "MEDIAN_TARGET_COVERAGE": {"condition": None},
                "PCT_TARGET_BASES_50X": {"condition": None},
                "PCT_TARGET_BASES_100X": {"condition": None},
                "PCT_TARGET_BASES_250X": {"condition": None},
                "PCT_TARGET_BASES_500X": {"condition": None},
                "PCT_TARGET_BASES_1000X": {"condition": None},
                "FOLD_80_BASE_PENALTY": {"condition": {"norm": "lt", "threshold": 5.0}},
                "PCT_OFF_BAIT": {"condition": {"norm": "lt", "threshold": 0.45}},
            },
        },
        "wgs": {
            "multiqc_picard_insertSize.json": {"MEAN_INSERT_SIZE": {"condition": None}},
            "multiqc_picard_dups.json": {"PERCENT_DUPLICATION": {"condition": None}},
        },
    }
}

METRICS_TO_DELIVER = {
    "targeted": [
        "MEAN_INSERT_SIZE",
        "PERCENT_DUPLICATION",
        "MEAN_TARGET_COVERAGE",
        "MEDIAN_TARGET_COVERAGE",
        "FOLD_80_BASE_PENALTY",
        "PCT_OFF_BAIT",
    ],
    "wgs": [
        "MEAN_INSERT_SIZE",
        "PERCENT_DUPLICATION",
        "FOLD_80_BASE_PENALTY",
    ],
}
