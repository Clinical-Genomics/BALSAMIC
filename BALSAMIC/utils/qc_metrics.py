from BALSAMIC.utils.constants import METRICS
from BALSAMIC.utils.models import QCCheckModel


def get_qc_metrics(analysis_path, sequencing_type):
    """Reads picard metrics for a particular case and returns a new dict

    Args:
        analysis_path: analysis result path e.g. /path/case_name/analysis
        sequencing_type : NGS approach
            targeted : if capture kit was used to enrich specific genomic regions
            wgs : if whole genome sequencing was performed

    Returns:
        qc_json: quality control summarised metrics

    """

    qc_check_model = QCCheckModel(
        analysis_path=analysis_path,
        sequencing_type=sequencing_type,
        qc_metrics=METRICS['qc']
    )

    qc_json = qc_check_model.get_metrics_json

    return qc_json
