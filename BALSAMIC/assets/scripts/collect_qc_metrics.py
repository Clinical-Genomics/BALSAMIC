#!/usr/bin/env python
import json
import os

import click
import yaml

from BALSAMIC.constants.quality_check_reporting import METRICS

from BALSAMIC.utils.models import MetricModel


@click.command(
    short_help="Extract the manually specified QC metrics",
)
@click.argument("output_path", type=click.Path(exists=False), required=True)
@click.argument("multiqc_data_path", type=click.Path(exists=True), required=True)
@click.argument("sequencing_type", required=True)
@click.argument("capture_kit", required=True)
def collect_qc_metrics(output_path, multiqc_data_path, sequencing_type, capture_kit):
    """Extracts the requested metrics and saves them to the output path"""

    with open(output_path, "w") as fn:
        yaml.dump(
            get_multiqc_metrics(multiqc_data_path, sequencing_type, capture_kit),
            fn,
            sort_keys=False,
            default_flow_style=False,
        )


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


def get_qc_available_panel_beds(metrics):
    """Returns available panel bed file names"""
    available_beds = []

    for k in metrics:
        if k != "default":
            available_beds.append(k)

    return available_beds


def get_requested_metrics(metrics, sequencing_type, capture_kit):
    """Parses the requested metrics and returns them as a dictionary"""

    requested_metrics = metrics[sequencing_type]
    if capture_kit:
        requested_metrics = metrics[sequencing_type]["default"]
        if capture_kit in get_qc_available_panel_beds(metrics[sequencing_type]):
            requested_metrics.update(metrics[sequencing_type][capture_kit])

    return requested_metrics


def get_multiqc_metrics(multiqc_data_path, sequencing_type, capture_kit):
    """Extracts the requested metrics from a multiqc JSON file"""

    with open(multiqc_data_path, "r") as f:
        multiqc_data = json.load(f)

    requested_metrics = get_requested_metrics(METRICS, sequencing_type, capture_kit)

    def extract(data, output_metrics, sample=None, source=None):
        """Recursively fetch metrics data from a nested multiqc JSON"""

        if isinstance(data, dict):
            for k in data:
                if "umi" not in k:
                    if k in requested_metrics:
                        output_metrics.append(
                            MetricModel(
                                id=sample.split("_")[1],
                                input=get_multiqc_data_source(
                                    multiqc_data, sample, source
                                ),
                                name=k,
                                step=source,
                                value=data[k],
                                condition=requested_metrics[k]["condition"],
                            ).dict()
                        )
                    extract(data[k], output_metrics, k, sample)

        return output_metrics

    return extract(multiqc_data["report_saved_raw_data"], [])


if __name__ == "__main__":
    collect_qc_metrics()
