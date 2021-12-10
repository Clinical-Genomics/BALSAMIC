import yaml

from BALSAMIC.utils.models import MetricValidationModel


def get_qc_metrics(yaml_path):
    """Retrieves the required and extracted QC metrics from a yaml file"""

    with open(yaml_path, "r") as fn:
        requested_metrics = yaml.load(fn, Loader=yaml.SafeLoader)

    return requested_metrics


def get_qc_metric_value(metrics, sample_id, metric_name):
    """Extracts the metrics value associated to a specific sample_id and metric_name"""

    for metric in metrics:
        if metric["id"] == sample_id and metric["name"] == metric_name:
            return metric["value"]

    return None


def validate_qc_metrics(metrics):
    """Returns a set of validated QC metrics"""

    return MetricValidationModel(metrics=metrics).dict()["metrics"]
