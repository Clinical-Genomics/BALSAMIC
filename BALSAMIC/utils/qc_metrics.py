import json
import os

from BALSAMIC.constants.quality_check_reporting import (
    METRICS,
    METRICS_TO_DELIVER,
)
from BALSAMIC.utils.models import QCCheckModel, DeliveryMetricModel


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


def get_qc_metrics_json(analysis_path, sequencing_type):
    """Extracts the metrics of interest and returns them as a json object"""
    qc_check_model = QCCheckModel.parse_obj(
        {"metrics": get_qc_metrics_dict(analysis_path, METRICS["qc"][sequencing_type])}
    )

    return qc_check_model.get_json


def get_multiqc_data_source(data, sample, source_name):
    """Extracts the metrics data source associated with sample and source names"""

    # Splits multiqc_picard_dups into ['multiqc', 'picard', 'dup'] in order to retrieve the
    # ["report_data_sources"]["Picard"]["DuplicationMetrics"] values from multiqc_data.json
    source = source_name[:-1].split("_")

    # Nested json fetching
    for source_tool in data["report_data_sources"]:
        for source_step in data["report_data_sources"][source_tool]:
            if (
                source[1].lower() in source_tool.lower()
                and source[2].lower() in source_step.lower()
            ):
                try:
                    return os.path.basename(
                        data["report_data_sources"][source_tool][source_step][sample]
                    )
                except KeyError:
                    # Deletes par orientation information from the sample name (insertSize metrics)
                    sample = sample.rsplit("_", 1)[0]

                    return os.path.basename(
                        data["report_data_sources"][source_tool][source_step][sample]
                    )


def extract_metrics_for_delivery(analysis_path, sequencing_type):
    """Extracts the output metrics to be delivered"""
    with open(
        os.path.join(analysis_path, "qc", "multiqc_data", "multiqc_data.json"), "r"
    ) as f:
        raw_data = json.load(f)

    def extract(data, output_metrics, sample=None, source=None):
        """Recursively fetch metrics information from nested multiQC JSON"""
        if isinstance(data, dict):
            for k in data:
                if "umi" not in k:
                    if k in METRICS_TO_DELIVER[sequencing_type]:
                        output_metrics.append(
                            DeliveryMetricModel(
                                id=sample.split("_")[1],
                                input=get_multiqc_data_source(raw_data, sample, source),
                                name=k,
                                step=source,
                                value=data[k],
                            ).dict()
                        )
                    extract(data[k], output_metrics, k, sample)

        return output_metrics

    return extract(raw_data["report_saved_raw_data"], [])
