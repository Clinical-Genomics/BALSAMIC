#!/usr/bin/env python
import json
import os
from pathlib import Path
from typing import List, Union

import click
import yaml

from BALSAMIC.constants.qc_metrics import METRICS
from BALSAMIC.utils.models import MetricModel


@click.command(
    short_help="Extract the manually specified QC metrics",
)
@click.argument("output_path", type=click.Path(exists=False), required=True)
@click.argument("multiqc_data_path", type=click.Path(exists=True), required=True)
@click.argument("counts_path", nargs=-1, type=click.Path(exists=True), required=False)
@click.argument("sequencing_type", required=True)
@click.argument("capture_kit", required=True)
def collect_qc_metrics(
    output_path: Path,
    multiqc_data_path: Path,
    counts_path: List[Path],
    sequencing_type: str,
    capture_kit: str,
):
    """Extracts the requested metrics from a JSON multiqc file and saves them to a YAML file

    Args:
        output_path: Path; destination path for the extracted YAML formatted metrics
        multiqc_data_path: Path; multiqc JSON path from which the metrics will be extracted
        counts_path: Path; list of variant caller specific files containing the number of variants
        sequencing_type: str; analysis sequencing type
        capture_kit: str; capture kit used for targeted analysis ("None" for WGS)
    """

    # MultiQC metrics
    metrics = get_multiqc_metrics(
        multiqc_data_path, sequencing_type, capture_kit_resolve_type(capture_kit)
    )

    # Number of variants
    for count in counts_path:
        metrics += get_variant_metrics(count)

    with open(output_path, "w") as fn:
        yaml.dump(
            metrics,
            fn,
            sort_keys=False,
            default_flow_style=False,
        )


def capture_kit_resolve_type(capture_kit: str):
    """Resolves the capture_kit type (NoneType or String)"""

    if capture_kit == "None":
        return None

    return capture_kit


def get_multiqc_data_source(multiqc_data: dict, sample: str, tool: str) -> str:
    """Extracts the metrics data source associated with a specific sample and tool

    Args:
        multiqc_data: dict; raw data from the multiqc_data.json file
        sample: str; sample ID
        tool: str; QC analysis tools applied during the workflow (e.g. "multiqc_picard_dups")

    Returns:
        A source file that was used to produce a specific metric
    """

    if tool == "multiqc_general_stats":
        subtool_name = ["multiqc", "FastQC", "all_sections"]
    else:
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


def get_qc_supported_capture_kit(capture_kit, metrics: List[str]) -> str:
    """Returns a BALSAMIC supported panel bed name associated to a specific capture_kit parameter"""
    available_panel_beds = []

    for k in metrics:
        if k != "default":
            available_panel_beds.append(k)

    return next((i for i in available_panel_beds if i in capture_kit), None)


def get_requested_metrics(
    metrics: dict, analysis_type: str, capture_kit: Union[str, None]
) -> dict:
    """Parses the defined and requested metrics and returns them as a dictionary"""

    requested_metrics = metrics[analysis_type]
    if capture_kit:
        requested_metrics = metrics[analysis_type]["default"]
        supported_capture_kit = get_qc_supported_capture_kit(
            capture_kit, metrics[analysis_type]
        )
        if supported_capture_kit:
            requested_metrics.update(metrics[analysis_type][supported_capture_kit])

    return requested_metrics


def get_multiqc_metrics(
    multiqc_data_path: Path, sequencing_type: str, capture_kit: Union[str, None]
) -> list:
    """Extracts and returns the requested metrics from a multiqc JSON file"""

    with open(multiqc_data_path, "r") as f:
        multiqc_data = json.load(f)

    requested_metrics = get_requested_metrics(METRICS, sequencing_type, capture_kit)

    def extract(data, output_metrics, sample=None, source=None):
        """Recursively fetch metrics data from a nested multiqc JSON"""

        if isinstance(data, dict):
            for k in data:
                # Ignore UMI and reverse reads metrics
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


def extract_number_variants(counts: list) -> dict:
    """Formats the number of SNPs, Indels, and total number of sites"""

    variant_metrics = dict()

    for count in counts:
        count = count.split(":")
        if len(count) > 1:
            variant_metrics.update(
                {count[0].strip().upper().replace(" ", "_"): int(count[1].strip())}
            )

    return variant_metrics


def get_variant_metrics(counts_path: list) -> list:
    """Retrieves the variant metrics and returns them as a MetricModel list"""

    output_metrics = list()

    with open(counts_path, "r") as input_file:
        counts = input_file.read().split("\n")

    variant_metrics = extract_number_variants(counts)
    requested_metrics = get_requested_metrics(METRICS, "variants", None)
    for metric in requested_metrics:
        output_metrics.append(
            MetricModel(
                id=os.path.basename(counts_path).split(".")[2],  # case_id
                input=os.path.basename(counts_path),
                name=metric,
                step="collect_custom_qc_metrics",
                value=variant_metrics[metric],
                condition=requested_metrics[metric]["condition"],
            ).dict()
        )

    return output_metrics


if __name__ == "__main__":
    collect_qc_metrics()
