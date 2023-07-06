"""QC metrics utility methods."""
from BALSAMIC.models.metrics import MetricValidation


def validate_qc_metrics(metrics: dict) -> dict:
    """Returns a set of validated QC metrics."""
    return MetricValidation(metrics=metrics).dict()["metrics"]
