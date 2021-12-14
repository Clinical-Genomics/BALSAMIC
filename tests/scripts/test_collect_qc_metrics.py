import copy
import json
import os
from pathlib import Path

from BALSAMIC.assets.scripts.collect_qc_metrics import (
    get_multiqc_data_source,
    get_multiqc_metrics,
    collect_qc_metrics,
    get_qc_available_panel_beds,
    merge_dicts,
    get_requested_metrics,
)


def test_get_qc_available_panel_beds(qc_requested_metrics):
    """test extraction of capture kits available for analysis"""

    # GIVEN an expected output
    expected_output = ["panel_1.bed", "panel_2.bed"]

    # WHEN calling the function
    available_panel_beds = get_qc_available_panel_beds(qc_requested_metrics["targeted"])

    # THEN check if the extracted bed file names correspond to the expected ones
    assert available_panel_beds == expected_output


def test_merge_dicts(qc_requested_metrics):
    """test dictionary merging and requirements overwriting by panel BED specific conditions"""

    # GIVEN an expected output
    expected_output = {
        "METRIC_1": {"condition": {"norm": "gt", "threshold": 1}},
        "METRIC_2": {"condition": {"norm": "gt", "threshold": 2}},
        "METRIC_4": {"condition": {"norm": "gt", "threshold": 4}},
    }

    # WHEN calling the function
    merged_dict = merge_dicts(
        copy.deepcopy(qc_requested_metrics["targeted"]["default"]),
        copy.deepcopy(qc_requested_metrics["targeted"]["panel_2.bed"]),
    )

    # THEN check if the extracted output corresponds to the merged dictionary
    assert merged_dict.items() == expected_output.items()


def test_get_requested_metrics_targeted(qc_requested_metrics):
    """test retrieval of the requested targeted metrics"""

    # GIVEN a sequencing type and a capture kit
    seq_type = "targeted"
    capture_kit = "panel_1.bed"

    # GIVEN the expected output
    expected_output = {
        "METRIC_1": {"condition": None},
        "METRIC_2": {"condition": {"norm": "gt", "threshold": 2}},
        "METRIC_3": {"condition": {"norm": "gt", "threshold": 3}},
    }

    # WHEN calling the function
    requested_metrics = get_requested_metrics(
        qc_requested_metrics, seq_type, capture_kit
    )

    # THEN check if the requested targeted metrics are correctly retrieved
    assert requested_metrics.items() == expected_output.items()


def test_get_requested_metrics_wgs(qc_requested_metrics):
    """test extraction of the requested WGS metrics"""

    # GIVEN a sequencing type and a capture kit
    seq_type = "wgs"
    capture_kit = None

    # GIVEN the expected output
    expected_output = {
        "METRIC_1": {"condition": {"norm": "gt", "threshold": 1}},
    }

    # WHEN calling the function
    requested_metrics = get_requested_metrics(
        qc_requested_metrics, seq_type, capture_kit
    )

    # THEN check if the requested metrics are WGS specific
    assert requested_metrics.items() == expected_output.items()


def test_get_multiqc_data_source(analysis_path):
    """test multiqc source extraction from multiqc_data.json analysis file"""

    # GIVEN input parameters and the multiqc data
    sample = "concatenated_tumor_XXXXXX_R"
    source_name_hs_metrics = "multiqc_picard_HsMetrics"
    source_name_dup = "multiqc_picard_dups"

    with open(
        os.path.join(analysis_path, "qc", "multiqc_data", "multiqc_data.json"), "r"
    ) as f:
        multiqc_data = json.load(f)

    # GIVEN an expected output
    source_hs_metrics = "concatenated_tumor_XXXXXX_R.sorted.mrkdup.hsmetric"
    source_dup = "concatenated_tumor_XXXXXX_R.sorted.mrkdup.txt"

    # WHEN extracting the source of a specific sample and collection of metrics
    out_source_hs_metrics = get_multiqc_data_source(
        multiqc_data, sample, source_name_hs_metrics
    )
    out_source_dup = get_multiqc_data_source(multiqc_data, sample, source_name_dup)

    # THEN check if the extracted source names correspond to the expected ones
    assert source_hs_metrics == out_source_hs_metrics
    assert source_dup == out_source_dup


def test_get_multiqc_metrics(analysis_path):
    """test metrics retrieval from the multiqc_data.json file"""

    # GIVEN a sequencing type and a capture kit
    seq_type = "targeted"
    capture_kit = "lymphoma_6.1_hg19_design.bed"

    # GIVEN an expected output
    n_metrics = 11  # Number of expected metric

    hs_metric = {
        "header": None,
        "id": "tumor",
        "input": "concatenated_tumor_XXXXXX_R.sorted.mrkdup.hsmetric",
        "name": "MEDIAN_TARGET_COVERAGE",
        "step": "multiqc_picard_HsMetrics",
        "value": 2393.0,
        "condition": {"norm": "gt", "threshold": 1000.0},
    }

    ins_size_metric = {
        "header": None,
        "id": "tumor",
        "input": "concatenated_tumor_XXXXXX_R.sorted.insertsizemetric",
        "name": "MEAN_INSERT_SIZE",
        "step": "multiqc_picard_insertSize",
        "value": 201.813054,
        "condition": None,
    }

    dups_metric = {
        "header": None,
        "id": "tumor",
        "input": "concatenated_tumor_XXXXXX_R.sorted.mrkdup.txt",
        "name": "PERCENT_DUPLICATION",
        "step": "multiqc_picard_dups",
        "value": 0.391429,
        "condition": None,
    }

    # WHEN calling the function
    metrics = get_multiqc_metrics(
        os.path.join(analysis_path, "qc", "multiqc_data", "multiqc_data.json"),
        seq_type,
        capture_kit,
    )

    # THEN check if the metrics are correctly retrieved
    assert len(metrics) == n_metrics
    assert hs_metric in metrics
    assert ins_size_metric in metrics
    assert dups_metric in metrics


def test_get_multiqc_metrics_filtering_umi(analysis_path):
    """tests that UMI data is filtered out when extracting metrics"""

    # GIVEN a sequencing type and a capture kit
    seq_type = "targeted"
    capture_kit = "lymphoma_6.1_hg19_design.bed"

    # WHEN calling the function
    metrics = get_multiqc_metrics(
        os.path.join(analysis_path, "qc", "multiqc_data", "multiqc_data.json"),
        seq_type,
        capture_kit,
    )

    # THEN check if the UMI samples are filtered out
    for metric in metrics:
        assert "umi" not in metric["input"]


def test_collect_qc_metrics(tmp_path, cli_runner):
    """tests qc metrics yaml file generation"""

    # GIVEN the output and multiqc metrics paths
    output_path = tmp_path / "tumor_metrics_deliverables.yaml"
    multiqc_data_path = (
        "tests/test_data/qc_files/analysis/qc/multiqc_data/multiqc_data.json"
    )

    # GIVEN a sequencing type and a capture kit
    seq_type = "targeted"
    capture_kit = "lymphoma_6.1_hg19_design.bed"

    # WHEN invoking the python script
    result = cli_runner.invoke(
        collect_qc_metrics,
        [str(output_path), multiqc_data_path, seq_type, capture_kit],
    )

    print(str(output_path))
    # THEN check if the YAML is correctly created and there are no errors
    assert result.exit_code == 0
    assert Path(output_path).exists()
