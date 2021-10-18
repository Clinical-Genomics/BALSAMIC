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

METRIC_FILES = {
    "picard_insertSize": "multiqc_picard_insertSize.json",
    "picard_dups": "multiqc_picard_dups.json",
    "picard_HsMetrics": "multiqc_picard_HsMetrics.json",
}

METRICS = {
    "qc": {
        "targeted": {
            "default": {
                METRIC_FILES["picard_insertSize"]: {
                    "MEAN_INSERT_SIZE": {"condition": None},
                },
                METRIC_FILES["picard_dups"]: {
                    "PERCENT_DUPLICATION": {"condition": None}
                },
                METRIC_FILES["picard_HsMetrics"]: {
                    "MEDIAN_TARGET_COVERAGE": {"condition": None},
                    "PCT_TARGET_BASES_50X": {"condition": None},
                    "PCT_TARGET_BASES_100X": {"condition": None},
                    "PCT_TARGET_BASES_250X": {"condition": None},
                    "PCT_TARGET_BASES_500X": {"condition": None},
                    "PCT_TARGET_BASES_1000X": {"condition": None},
                    "MEAN_TARGET_COVERAGE": {"condition": None},
                    "FOLD_80_BASE_PENALTY": {"condition": None},
                    "PCT_OFF_BAIT": {"condition": None},
                },
            },
            "gicfdna_3.1_hg19_design.bed": {
                METRIC_FILES["picard_HsMetrics"]: {
                    "MEAN_TARGET_COVERAGE": {
                        "condition": {"norm": "gt", "threshold": 500.0}
                    },
                    "FOLD_80_BASE_PENALTY": {
                        "condition": {"norm": "lt", "threshold": 1.5}
                    },
                    "PCT_OFF_BAIT": {"condition": {"norm": "lt", "threshold": 0.35}},
                }
            },
            "gmcksolid_4.1_hg19_design.bed": {
                METRIC_FILES["picard_HsMetrics"]: {
                    "MEAN_TARGET_COVERAGE": {
                        "condition": {"norm": "gt", "threshold": 500}
                    },
                    "FOLD_80_BASE_PENALTY": {
                        "condition": {"norm": "lt", "threshold": 1.7}
                    },
                    "PCT_OFF_BAIT": {"condition": {"norm": "lt", "threshold": 0.3}},
                }
            },
            "gmsmyeloid_5.2_hg19_design.bed": {
                METRIC_FILES["picard_HsMetrics"]: {
                    "MEAN_TARGET_COVERAGE": {
                        "condition": {"norm": "gt", "threshold": 1000}
                    },
                    "FOLD_80_BASE_PENALTY": {
                        "condition": {"norm": "lt", "threshold": 1.5}
                    },
                    "PCT_OFF_BAIT": {"condition": {"norm": "lt", "threshold": 0.4}},
                }
            },
            "lymphoma_6.1_hg19_design.bed": {
                METRIC_FILES["picard_HsMetrics"]: {
                    "MEAN_TARGET_COVERAGE": {
                        "condition": {"norm": "gt", "threshold": 1000}
                    },
                    "FOLD_80_BASE_PENALTY": {
                        "condition": {"norm": "lt", "threshold": 1.5}
                    },
                    "PCT_OFF_BAIT": {"condition": {"norm": "lt", "threshold": 0.35}},
                }
            },
            "gmslymphoid_7.1_hg19_design.bed": {
                METRIC_FILES["picard_HsMetrics"]: {
                    "MEAN_TARGET_COVERAGE": {
                        "condition": {"norm": "gt", "threshold": 1000}
                    },
                    "FOLD_80_BASE_PENALTY": {
                        "condition": {"norm": "lt", "threshold": 1.5}
                    },
                    "PCT_OFF_BAIT": {"condition": {"norm": "lt", "threshold": 0.35}},
                }
            },
            "twistexomerefseq_9.1_hg19_design.bed": {
                METRIC_FILES["picard_HsMetrics"]: {
                    "MEAN_TARGET_COVERAGE": {
                        "condition": {"norm": "gt", "threshold": 100}
                    },
                    "FOLD_80_BASE_PENALTY": {
                        "condition": {"norm": "lt", "threshold": 1.8}
                    },
                    "PCT_OFF_BAIT": {"condition": {"norm": "lt", "threshold": 0.25}},
                }
            },
        },
        "wgs": {
            METRIC_FILES["picard_insertSize"]: {
                "MEAN_INSERT_SIZE": {"condition": None}
            },
            METRIC_FILES["picard_dups"]: {"PERCENT_DUPLICATION": {"condition": None}},
        },
    }
}
