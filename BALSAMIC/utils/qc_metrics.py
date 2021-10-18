import json
import os

from BALSAMIC.constants.quality_check_reporting import METRICS
from BALSAMIC.utils.models import QCCheckModel


def get_qc_available_panel_beds(metrics):
    """Returns available panel beds file names for QC validation"""
    available_beds = []

    for k in metrics:
        if k != "default":
            available_beds.append(k)

    return available_beds


def merge_dicts(*dicts):
    """Merges multiple dictionaries integrating by common keys"""
    merged_dict = {}

    for d in dicts:
        for key in d:
            try:
                # Overwrites the default values with panel specific ones
                merged_dict[key].update(d[key])
            except KeyError:
                merged_dict[key] = d[key]

    return merged_dict


def read_metrics(analysis_path, file_name):
    """Extracts all the metrics from a specific QC file"""
    with open(os.path.join(analysis_path, "qc", "multiqc_data", file_name), "r") as f:
        raw_metrics = json.load(f)

    # Ignore the metrics associated with UMIs
    filtered_raw_metrics = {
        sample_name: metrics
        for sample_name, metrics in raw_metrics.items()
        if "umi" not in sample_name
    }

    return filtered_raw_metrics


def update_metrics_dict(sample_id, metric, value, metrics_dict):
    """Appends a {metric, value, condition} object to a dictionary"""
    sample_name = "_".join([sample_id.split("_")[0], sample_id.split("_")[1]])

    if sample_name not in metrics_dict:
        metrics_dict[sample_name] = []

    metrics_dict[sample_name].append(
        {"name": metric[0], "value": value, "condition": metric[1]["condition"]}
    )

    return metrics_dict


def get_qc_metrics_dict(analysis_path, requested_metrics):
    """Returns a dictionary of the requested QC metrics along with their values and filtering conditions"""
    metrics_dict = {}

    # Loop through MultiQC json files
    for file_name, metrics in requested_metrics.items():
        raw_metrics = read_metrics(analysis_path, file_name)
        for j in raw_metrics:
            for k in metrics.items():
                metrics_dict = update_metrics_dict(
                    j, k, raw_metrics[j][k[0]], metrics_dict
                )
    return metrics_dict


def get_qc_metrics_json(analysis_path, sequencing_type, panel_bed):
    """Extracts the metrics of interest and returns them as a json object"""
    if sequencing_type != "wgs" and panel_bed in get_qc_available_panel_beds(
        METRICS["qc"][sequencing_type]
    ):
        metrics = merge_dicts(
            METRICS["qc"][sequencing_type]["default"],
            METRICS["qc"][sequencing_type][panel_bed],
        )
    elif sequencing_type != "wgs":
        metrics = METRICS["qc"][sequencing_type]["default"]
    else:
        metrics = METRICS["qc"][sequencing_type]

    qc_check_model = QCCheckModel.parse_obj(
        {"metrics": get_qc_metrics_dict(analysis_path, metrics)}
    )

    return qc_check_model.get_json
