from BALSAMIC.utils.models import MetricValidationModel


def validate_qc_metrics(metrics: dict) -> dict:
    """Returns a set of validated QC metrics"""

    return MetricValidationModel(metrics=metrics).dict()["metrics"]
