"""QC metrics utility methods."""
from BALSAMIC.models.metrics import MetricValidation


def validate_qc_metrics(metrics: dict, LOG) -> dict:
    """Returns a set of validated QC metrics."""
    return MetricValidation(metrics=metrics).model_dump()["metrics"]
