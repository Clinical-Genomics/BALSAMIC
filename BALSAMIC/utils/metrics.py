"""QC metrics utility methods."""
from BALSAMIC.models.metrics_model import MetricValidationModel


def validate_qc_metrics(metrics: dict) -> dict:
    """Returns a set of validated QC metrics."""
    return MetricValidationModel(metrics=metrics).dict()["metrics"]
