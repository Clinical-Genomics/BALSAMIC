import json
import os

from BALSAMIC.utils.constants import METRICS
from BALSAMIC.utils.models import QCExtractionModel

MULTIPLE_JSON = [
    {
        "file_name": "multiqc_picard_insertSize.json",
        "required_metrics": ["MEAN_INSERT_SIZE"],
    },
    {
        "file_name": "multiqc_picard_dups.json",
        "required_metrics": ["PERCENT_DUPLICATION"],
    },
    {
        "file_name": "multiqc_picard_HsMetrics.json",
        "required_metrics": [
            "MEAN_TARGET_COVERAGE",
            "MEDIAN_TARGET_COVERAGE",
            "PCT_TARGET_BASES_50X",
            "PCT_TARGET_BASES_100X",
            "PCT_TARGET_BASES_250X",
            "PCT_TARGET_BASES_500X",
            "PCT_TARGET_BASES_1000X",
            "FOLD_80_BASE_PENALTY",
        ],
    },
]


def get_qc_metrics(analysis_path):
    """Reads picard metrics for a particular case and returns a new dict

    Args:
        analysis_path: analysis result path e.g. /path/case_name/analysis

    Returns:
        qc_data: a dictionary of extracted metrics

    """
    qc_data = {}

    # Loop through json files
    for single_json in MULTIPLE_JSON:
        file_name = os.path.join(
            analysis_path, "qc", "multiqc_data", single_json["file_name"]
        )
        with open(file_name, "r") as f:
            json_file = json.load(f)
        for k in single_json["required_metrics"]:
            for j in json_file.keys():
                if "umi" not in j:
                    sampleid = j.split("_")[1]
                    if sampleid not in qc_data:
                        qc_data[sampleid] = {}
                    qc_data[sampleid][k] = json_file[j][k]

    return qc_data


def get_qc_metrics_json(analysis_path, sequencing_type):
    """Extracts the metrics of interest and returns them as a json object

    Args:
        analysis_path: analysis result path e.g. /path/case_name/analysis
        sequencing_type : NGS approach
            targeted : if capture kit was used to enrich specific genomic regions
            wgs : if whole genome sequencing was performed

    Returns:
        qc_metrics_json: quality control summarised metrics

    """

    qc_extraction_model = QCExtractionModel(
        analysis_path=analysis_path,
        sequencing_type=sequencing_type,
        qc_attributes=METRICS["qc"],
    )

    qc_metrics_json = json.dumps(
        qc_extraction_model.get_metrics, indent=4, sort_keys=True
    )

    return qc_metrics_json
