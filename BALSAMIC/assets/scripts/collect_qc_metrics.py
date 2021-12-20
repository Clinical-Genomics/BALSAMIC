#!/usr/bin/env python
import json
import os
from pathlib import Path
from typing import List, Union

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
def collect_qc_metrics(
    output_path: Path,
    multiqc_data_path: Path,
    sequencing_type: str,
    capture_kit: Union[str, None],
):
    """Extracts the requested metrics from a JSON multiqc file and saves them to a YAML file

    Args:
        output_path: Path; destination path for the extracted YAML formatted metrics
        multiqc_data_path: Path; multiqc JSON path from which the metrics will be extracted
        sequencing_type: str; analysis sequencing type
        capture_kit: str; capture kit used for targeted analysis (None for WGS)
    """

    with open(output_path, "w") as fn:
        yaml.dump(
            get_multiqc_metrics(multiqc_data_path, sequencing_type, capture_kit),
            fn,
            sort_keys=False,
            default_flow_style=False,
        )


def get_multiqc_data_source(multiqc_data: dict, sample: str, tool: str) -> str:
    """Extracts the metrics data source associated with a specific sample and tool

    Args:
        multiqc_data: dict; raw data from the multiqc_data.json file
        sample: str; sample ID
        tool: str; QC analysis tools applied during the workflow (e.g. "multiqc_picard_dups")

    Returns:
        A source file that was used to produce a specific metric
    """

    # Use case: splits multiqc_picard_dups into ['multiqc', 'picard', 'dup'] in order to retrieve the
    # ["report_data_sources"]["Picard"]["DuplicationMetrics"] values from multiqc_data.json
    subtool_name = tool[:-1].split("_")

    # Nested json fetching
    for source_tool in multiqc_data["report_data_sources"]:
        # source_tool: Picard, fastp, FastQC, etc.
        for source_subtool in multiqc_data["report_data_sources"][source_tool]:
            # source_subtool (for Picard): AlignmentSummaryMetrics, HsMetrics, DuplicationMetric, etc.
            if (
                subtool_name[1].lower() in source_tool.lower()
                and subtool_name[2].lower() in source_subtool.lower()
            ):
                try:
                    return os.path.basename(
                        multiqc_data["report_data_sources"][source_tool][
                            source_subtool
                        ][sample]
                    )
                except KeyError:
                    # Deletes par orientation information from the sample name (insertSize metrics)
                    sample = sample.rsplit("_", 1)[0]
                    return os.path.basename(
                        multiqc_data["report_data_sources"][source_tool][
                            source_subtool
                        ][sample]
                    )


def get_qc_available_panel_beds(metrics: List[str]) -> List[str]:
    """Returns available panel bed file names from a list of requested metrics"""
    available_beds = []

    for k in metrics:
        if k != "default":
            available_beds.append(k)

    return available_beds


def get_requested_metrics(
    metrics: dict, sequencing_type: str, capture_kit: Union[str, None]
) -> dict:
    """Parses the defined and requested metrics and returns them as a dictionary"""

    requested_metrics = metrics[sequencing_type]
    if capture_kit:
        requested_metrics = metrics[sequencing_type]["default"]
        if capture_kit in get_qc_available_panel_beds(metrics[sequencing_type]):
            requested_metrics.update(metrics[sequencing_type][capture_kit])

    return requested_metrics


def get_multiqc_metrics(
    multiqc_data_path: Path, sequencing_type: str, capture_kit: Union[str, None]
) -> dict:
    """Extracts the requested metrics from a multiqc JSON file and returns them as a dictionary"""

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
