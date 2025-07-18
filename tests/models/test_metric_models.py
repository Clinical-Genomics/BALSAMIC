"""Tests for the QC metrics related methods."""
import copy
from typing import Any, Dict
import pytest

from BALSAMIC.models.metrics import MetricValidation, Metric, MetricCondition
from BALSAMIC.assets.scripts.collect_qc_metrics import get_sex_check_metrics


def test_metric_condition():
    """Test MetricCondition attributes parsing."""

    # GIVEN input attributes
    metric_condition: Dict[str, Any] = {"norm": "gt", "threshold": 1}

    # WHEN building the metric condition model
    metric_model: MetricCondition = MetricCondition(**metric_condition)

    # THEN assert retrieved values from the created model
    assert metric_model.model_dump() == metric_condition


def test_metric_pass_validation():
    """Test Metric attributes parsing."""

    # GIVEN input attributes
    metrics: Dict[str, Any] = {
        "header": None,
        "id": "ACC1",
        "input": "ACC1.sorted.mrkdup.hsmetric",
        "name": "MEDIAN_TARGET_COVERAGE",
        "step": "multiqc_picard_HsMetrics",
        "value": 2393.0,
        "condition": {"norm": "gt", "threshold": 1000.0},
    }

    # WHEN building the metric model
    metric_model: Metric = Metric(**metrics)

    # THEN assert retrieved values from the created model
    assert metric_model.model_dump() == metrics


def test_metric_fail_validation():
    """Test Metric behaviour for an incorrect input."""

    # GIVEN an invalid input
    invalid_input: Dict[str, Any] = {"header": None, "id": "ACC1"}

    # THEN the model raises an error due to an incomplete input
    with pytest.raises(ValueError) as input_exc:
        Metric(**invalid_input)
    assert "Field required" in str(input_exc.value)


def test_metric_validation_pass(qc_extracted_metrics: dict):
    """Test MetricValidation attribute parsing and positive validation."""

    # WHEN building the MetricValidation model
    model: MetricValidation = MetricValidation(metrics=qc_extracted_metrics)

    # THEN assert retrieved values from the created model
    assert model.model_dump()["metrics"] == qc_extracted_metrics


def test_metric_validation_fail(qc_extracted_metrics: dict):
    """Test MetricValidation for an overly restrictive metric condition."""

    # GIVEN input attributes with a value that does not meet the filtering condition
    metrics: dict = copy.deepcopy(qc_extracted_metrics)
    metrics[4]["value"] = 2.0  # GC_DROPOUT set to 2.0 (failing condition)

    # THEN check that the model filters the metric according to its norm
    with pytest.raises(ValueError) as val_exc:
        MetricValidation(metrics=metrics)
    assert (
        f"QC metric {metrics[4]['name']}: {metrics[4]['value']} validation has failed. "
        f"(Condition: {metrics[4]['condition']['norm']} {metrics[4]['condition']['threshold']}, ID: {metrics[4]['id']})"
        in str(val_exc.value)
    )


def test_multiple_metric_validation_fail(qc_extracted_metrics: dict):
    """Test MetricValidation for multiple metrics with failing conditions."""

    # GIVEN input attributes that does not meet the specified conditions
    metrics: dict = copy.deepcopy(qc_extracted_metrics)
    metrics[4]["value"] = 2.0  # GC_DROPOUT set to 2.0 (failing condition)
    metrics[8]["value"] = 0.5  # PCT_TARGET_BASES_500X set to 50% (failing condition)

    # THEN check that the model filters the metrics according to its norm
    with pytest.raises(ValueError) as val_exc:
        MetricValidation(metrics=metrics)
    assert "2 validation errors for MetricValidation" in str(val_exc.value)
    assert metrics[4]["name"] in str(val_exc.value)
    assert metrics[8]["name"] in str(val_exc.value)


def test_sex_check_metric_is_only_warned(
    tga_female_sex_prediction: str,
    config_dict: dict,
):
    """MetricValidation should no longer raise for COMPARE_PREDICTED_TO_GIVEN_SEX."""
    # GIVEN a female-prediction JSON and a male value in the config
    sex_metrics = get_sex_check_metrics(tga_female_sex_prediction, config_dict)

    # WHEN validating the metrics
    # THEN it should **not** raise; if it does, the test will fail automatically.
    MetricValidation(metrics=sex_metrics)


def test_validate_metric_type_error():
    """Test MetricValidation TypeError for incompatible variable-type comparisons."""

    # Create a Metric object where value and threshold are incompatible types
    metrics = [
        {
            "name": "Test Metric",
            "input": "fake_input_file",
            "step": "fake_step",
            "value": "string-value",  # Incompatible with numerical comparisons
            "condition": {"norm": "lt", "threshold": 10},
            "id": "caseid",
        }
    ]

    # Validate that the TypeError handling raises the correct ValueError
    with pytest.raises(ValueError, match=r"Type mismatch in QC metric .*"):
        MetricValidation(metrics=metrics)


def test_metric_validation_norm_fail(qc_extracted_metrics: dict):
    """Test MetricValidation ValueError raising for an operator that it is not accepted."""

    # GIVEN a metric with an incorrect norm attribute
    metrics: dict = copy.deepcopy(qc_extracted_metrics)
    metrics[4]["condition"]["norm"] = "lower"

    # THEN model raises an error due to a non accepted norm
    try:
        MetricValidation(metrics=metrics)
    except ValueError as key_exc:
        assert metrics[4]["condition"]["norm"] in str(key_exc)
