import os

from BALSAMIC.utils.qc_metrics import (
    get_qc_metrics,
    validate_qc_metrics,
    get_qc_metric_value,
)


def test_get_qc_metrics(metrics_yaml_path):
    """test requested metric extraction from the saved YAML file"""

    # GIVEN a metrics yaml file
    metrics_yaml_path = os.path.join(metrics_yaml_path)

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
    requested_metrics = get_qc_metrics(metrics_yaml_path)

    # THEN check if the metrics are correctly retrieved
    assert len(requested_metrics) == n_metrics
    assert hs_metric in requested_metrics
    assert ins_size_metric in requested_metrics
    assert dups_metric in requested_metrics


def test_get_qc_metric_value(qc_extracted_metrics):
    """test QC metric value extraction"""

    # GIVEN the input parameters
    sample_id = "tumor"
    metric_name = "MEDIAN_TARGET_COVERAGE"

    # GIVEN the expected value
    expected_value = 2393.0

    # WHEN calling the function
    metric_value = get_qc_metric_value(qc_extracted_metrics, sample_id, metric_name)

    # THEN check if the retrieved value corresponds to the expected one
    assert metric_value == expected_value


def test_validate_qc_metrics(qc_extracted_metrics):
    """test QC metric validation"""

    # WHEN calling the function
    validated_metrics_pass = validate_qc_metrics(qc_extracted_metrics)

    # THEN check if the obtained metrics are correctly parsed and validated
    assert validated_metrics_pass == qc_extracted_metrics
