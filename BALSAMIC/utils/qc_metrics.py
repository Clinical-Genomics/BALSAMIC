import json

from BALSAMIC.constants.quality_check_reporting import METRICS
from BALSAMIC.utils.models import QCExtractionModel


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
