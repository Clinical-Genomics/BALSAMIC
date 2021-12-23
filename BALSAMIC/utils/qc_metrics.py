from typing import Union

from BALSAMIC.utils.models import MetricValidationModel


def get_qc_metric_value(
    metrics: dict, sample_id: str, metric_name: str
) -> Union[float, None]:
    """Extracts the metrics value associated to a specific sample_id and metric_name"""

    for metric in metrics:
        if metric["id"] == sample_id and metric["name"] == metric_name:
            return metric["value"]

    return None


def validate_qc_metrics(metrics: dict) -> dict:
    """Returns a set of validated QC metrics"""

    return MetricValidationModel(metrics=metrics).dict()["metrics"]
