import os

from BALSAMIC.utils.qc_metrics import (
    validate_qc_metrics,
    get_qc_metric_value,
)


def test_get_qc_metric_value(qc_extracted_metrics):
    """test QC metric value extraction"""

    # GIVEN the input parameters
    sample_id = "tumor"
    metric_name = "MEDIAN_TARGET_COVERAGE"

    # GIVEN an expected value
    expected_value = 2393.0

    # WHEN calling the function
    metric_value = get_qc_metric_value(qc_extracted_metrics, sample_id, metric_name)

    # THEN check if the retrieved value corresponds to the expected one
    assert metric_value == expected_value


def test_get_qc_metric_value_invalid_metric(qc_extracted_metrics):
    """test QC metric value extraction for an invalid metric name"""

    # GIVEN the input parameters
    sample_id = "tumor"
    metric_name = "NOT_A_METRIC"

    # WHEN calling the function
    metric_value = get_qc_metric_value(qc_extracted_metrics, sample_id, metric_name)

    # THEN check if the retrieved value is None
    assert metric_value is None


def test_validate_qc_metrics(qc_extracted_metrics):
    """test QC metric validation"""

    # WHEN calling the function
    validated_metrics_pass = validate_qc_metrics(qc_extracted_metrics)

    # THEN check if the obtained metrics are correctly parsed and validated
    assert validated_metrics_pass == qc_extracted_metrics
