#!/usr/bin/env python
import os
import re
from pathlib import Path
from typing import List, Optional

import click
import yaml

from BALSAMIC.constants.metrics import METRICS
from BALSAMIC.models.metrics import Metric
from BALSAMIC.utils.io import read_json
from BALSAMIC.utils.rule import (
    get_analysis_type,
    get_capture_kit,
    get_sample_type_from_sample_name,
    get_sequencing_type,
)


@click.command(
    short_help="Extract the manually specified QC metrics",
)
@click.argument("config_path", type=click.Path(exists=True), required=True)
@click.argument("output_path", type=click.Path(exists=False), required=True)
@click.argument("multiqc_data_path", type=click.Path(exists=True), required=True)
@click.argument("counts_path", nargs=-1, type=click.Path(exists=True), required=False)
@click.option(
    "--sex-prediction-path",
    type=click.Path(exists=True),
    required=False,
    help="Path to sex prediction json (optional for balsamic-qc)",
)
def collect_qc_metrics(
    config_path: Path,
    output_path: Path,
    multiqc_data_path: Path,
    sex_prediction_path: Path,
    counts_path: List[Path],
):
    """Extracts the requested metrics from a JSON multiqc file and saves them to a YAML file

    Args:
        config_path: Path; case config file path
        output_path: Path; destination path for the extracted YAML formatted metrics
        multiqc_data_path: Path; multiqc JSON path from which the metrics will be extracted
        sex_prediction_path: Path; sex prediction JSON path from which sex prediction info will be extracted
        counts_path: Path; list of variant caller specific files containing the number of variants
    """

    config = read_json(config_path)
    multiqc_data = read_json(multiqc_data_path)

    # MultiQC metrics
    metrics = get_multiqc_metrics(config, multiqc_data)

    # Number of variants
    for count in counts_path:
        metrics += get_variant_metrics(count)

    # Sex check
    if sex_prediction_path:
        metrics += get_sex_check_metrics(sex_prediction_path, config)

    # Relatedness
    analysis_type = get_analysis_type(config)
    if analysis_type == "paired" and "Somalier" in multiqc_data["report_data_sources"]:
        metrics += get_relatedness_metrics(multiqc_data)

    with open(output_path, "w") as fn:
        yaml.dump(
            metrics,
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
                    # Deletes pair orientation information from the sample name (insertSize metrics)
                    sample = sample.rsplit("_", 1)[0]
                    return os.path.basename(
                        multiqc_data["report_data_sources"][source_tool][
                            source_subtool
                        ][sample]
                    )


def get_sex_check_metrics(sex_prediction_path: str, config: dict) -> list:
    """Retrieves the sex check metrics and returns them as a Metric list."""
    metric = "compare_predicted_to_given_sex"
    case_id: str = config["analysis"]["case_id"]
    sex_prediction: dict = read_json(sex_prediction_path)

    given_sex: str = config["analysis"]["gender"]

    sex_check_metrics = []

    for sample_type in ["tumor", "normal"]:
        if sample_type in sex_prediction:
            predicted_sex = sex_prediction[sample_type]["predicted_sex"]
            sex_prediction_metrics = Metric(
                id=f"{case_id}_{sample_type}",
                input=os.path.basename(sex_prediction_path),
                name=metric.upper(),
                step="sex_check",
                value=predicted_sex,
                condition={"norm": "eq", "threshold": given_sex},
            ).model_dump()
            sex_check_metrics.append(sex_prediction_metrics)

    return sex_check_metrics


def get_relatedness_metrics(multiqc_data: dict) -> list:
    """Retrieves the relatedness metrics and returns them as a Metric list."""
    source_tool = "Somalier"
    metric = "relatedness"
    step = "multiqc_somalier"

    # Somalier doesn't input the sample name but {$case_id}_{[NORMAL|TUMOR]},
    # where [NORMAL|TUMOR] is the sample group in the bam file
    for sample in multiqc_data["report_data_sources"][source_tool]["all_sections"]:
        if "*" in sample:
            data_source = os.path.basename(
                multiqc_data["report_data_sources"][source_tool]["all_sections"][sample]
            )
            metric_value = multiqc_data["report_saved_raw_data"][step][sample][metric]
            case_id = re.sub(r"_NORMAL.*|_TUMOR.*", "", sample)

            output_metrics = Metric(
                id=case_id,
                input=data_source,
                name=metric.upper(),
                step=step,
                value=metric_value,
                condition=METRICS["paired"]["RELATEDNESS"]["condition"],
            ).model_dump()

            return [output_metrics]


def get_qc_supported_capture_kit(capture_kit, metrics: List[str]) -> str:
    """Returns a BALSAMIC supported panel bed name associated to a specific capture_kit parameter"""
    available_panel_beds = []

    for k in metrics:
        if k != "default":
            available_panel_beds.append(k)

    return next(
        (
            i
            for i in available_panel_beds
            if re.search(rf"{re.escape(i)}(?=_\d)", capture_kit)
        ),
        None,
    )


def get_requested_metrics(config: dict, metrics: dict) -> dict:
    """Parses the defined and requested metrics and returns them as a dictionary"""

    sequencing_type = get_sequencing_type(config)
    capture_kit = get_capture_kit(config)

    requested_metrics = metrics[sequencing_type]
    if capture_kit:
        requested_metrics = metrics[sequencing_type]["default"]
        supported_capture_kit = get_qc_supported_capture_kit(
            capture_kit, metrics[sequencing_type]
        )
        if supported_capture_kit:
            requested_metrics.update(metrics[sequencing_type][supported_capture_kit])

    return requested_metrics


def get_metric_condition(
    config: dict, requested_metrics: dict, sample: str, metric: str
) -> Optional[dict]:
    """Returns a condition associated to a sample and sequencing type"""

    sequencing_type = get_sequencing_type(config)
    sample_type = get_sample_type_from_sample_name(config, sample)
    if not sample_type:
        # Deletes pair orientation information from the sample name (insertSize metrics)
        sample_type = get_sample_type_from_sample_name(config, sample.split("_")[0])

    req_metrics = requested_metrics[metric]["condition"]
    if sequencing_type == "wgs" and (
        (metric == "PCT_60X" and sample_type == "normal")
        or (metric == "MEDIAN_COVERAGE" and sample_type == "tumor")
    ):
        req_metrics = None

    return req_metrics


def get_sample_id(multiqc_key: str) -> str:
    """Returns extracted sample ID from MultiQC data JSON key.

    Example of possible sample-formats below from "report_saved_raw_data":
    tumor.ACCXXXXXX
    tumor.ACCXXXXXX_FR
    ACCXXXXXX_align_sort_HMYLNDSXX_ACCXXXXXX_S165_L001

    Returns
        str: The extracted sample ID with the ACCXXXXXX format.
    """
    if "_align_sort_" in multiqc_key:
        return multiqc_key.split("_")[0]
    return multiqc_key.split(".")[1].split("_")[0]


def get_multiqc_metrics(config: dict, multiqc_data: dict) -> list:
    """Extracts and returns the requested metrics from a multiqc JSON file"""

    requested_metrics = get_requested_metrics(config, METRICS)

    def extract(data, output_metrics, multiqc_key=None, source=None):
        """Recursively fetch metrics data from a nested multiqc JSON"""

        if isinstance(data, dict):
            for k in data:
                # Ignore UMI and reverse reads metrics
                if "umi" not in k:
                    if k in requested_metrics:
                        output_metrics.append(
                            Metric(
                                id=get_sample_id(multiqc_key),
                                input=get_multiqc_data_source(
                                    multiqc_data, multiqc_key, source
                                ),
                                name=k,
                                step=source,
                                value=data[k],
                                condition=get_metric_condition(
                                    config,
                                    requested_metrics,
                                    get_sample_id(multiqc_key),
                                    k,
                                ),
                            ).model_dump()
                        )
                    extract(data[k], output_metrics, k, multiqc_key)

        return output_metrics

    return extract(multiqc_data["report_saved_raw_data"], [])


def extract_number_variants(counts: list) -> dict:
    """Formats the number of SNPs, Indels, and total number of sites"""

    variant_metrics = dict()

    for count in counts:
        # Transforms string "Number of sites:   125" into a key value object {"NUMBER_OF_SITES": 125}
        count = count.split(":")
        if len(count) > 1:
            variant_metrics.update(
                {count[0].strip().upper().replace(" ", "_"): int(count[1].strip())}
            )

    return variant_metrics


def get_variant_metrics(counts_path: list) -> list:
    """Retrieves the variant metrics and returns them as a Metric list"""

    output_metrics = list()

    with open(counts_path, "r") as input_file:
        counts = input_file.read().split("\n")

    variant_metrics = extract_number_variants(counts)
    requested_metrics = METRICS["variants"]
    for metric in requested_metrics:
        output_metrics.append(
            Metric(
                id=os.path.basename(counts_path).split(".")[2],  # case_id
                input=os.path.basename(counts_path),
                name=metric,
                step="collect_custom_qc_metrics",
                value=variant_metrics[metric],
                condition=requested_metrics[metric]["condition"],
            ).model_dump()
        )

    return output_metrics


if __name__ == "__main__":
    collect_qc_metrics()
