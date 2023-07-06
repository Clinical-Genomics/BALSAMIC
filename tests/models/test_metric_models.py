"""Tests for the QC metrics related methods."""
import copy

import pytest

from BALSAMIC.models.metrics import (
    MetricValidation,
    Metric,
    MetricCondition,
)


def test_metric_condition_model():
    """test MetricCondition attributes parsing"""

    # GIVEN input attributes
    metric_condition = {"norm": "gt", "threshold": 1}

    # WHEN building the metric condition model
    metrics_model = MetricCondition(**metric_condition)

    # THEN assert retrieved values from the created model
    assert metrics_model.dict().items() == metric_condition.items()


def test_metric_model_pass_validation():
    """test Metric attributes parsing"""

    # GIVEN input attributes
    metrics = {
        "header": None,
        "id": "ACC1",
        "input": "ACC1.sorted.mrkdup.hsmetric",
        "name": "MEDIAN_TARGET_COVERAGE",
        "step": "multiqc_picard_HsMetrics",
        "value": 2393.0,
        "condition": {"norm": "gt", "threshold": 1000.0},
    }

    # WHEN building the metric model
    metric_model = Metric(**metrics)

    # THEN assert retrieved values from the created model
    assert metric_model.dict().items() == metrics.items()


def test_metric_model_duplication_refactoring():
    """test Metric duplications param refactoring"""

    # GIVEN input attributes
    metrics = {
        "header": None,
        "id": "ACC1",
        "input": "ACC1_R_1_fastqc.zip",
        "name": "FastQC_mqc-generalstats-fastqc-percent_duplicates",
        "step": "multiqc_general_stats",
        "value": 21.517800000611373,
        "condition": None,
    }

    # WHEN building the metric model
    metric_model = Metric(**metrics)

    # THEN assert retrieved values from the created model
    assert metric_model.name == "PERCENT_DUPLICATION_R1"


def test_metric_model_fail_validation():
    """test Metric behaviour for an incorrect input"""

    # GIVEN an invalid input
    invalid_input = {"header": None, "id": "ACC1"}

    # THEN the model raises an error due to an incomplete input
    with pytest.raises(ValueError) as input_exc:
        Metric(**invalid_input)
    assert "field required" in str(input_exc.value)


def test_metric_validation_model_pass(qc_extracted_metrics):
    """test MetricValidation attribute parsing and positive validation"""

    # WHEN building the MetricValidation model
    model = MetricValidation(metrics=qc_extracted_metrics)

    # THEN assert retrieved values from the created model
    assert model.dict()["metrics"] == qc_extracted_metrics


def test_metric_validation_model_fail(qc_extracted_metrics):
    """test MetricValidation for an overly restrictive metric condition"""

    # GIVEN input attributes with a value that does not meet the filtering condition
    metrics = copy.deepcopy(qc_extracted_metrics)
    metrics[4]["value"] = 2.0  # GC_DROPOUT set to 2.0 (failing condition)

    # THEN check that the model filters the metric according to its norm
    with pytest.raises(ValueError) as val_exc:
        MetricValidation(metrics=metrics)
    assert (
        f"QC metric {metrics[4]['name']}: {metrics[4]['value']} validation has failed. "
        f"(Condition: {metrics[4]['condition']['norm']} {metrics[4]['condition']['threshold']}, ID: {metrics[4]['id']})"
        in str(val_exc.value)
    )


def test_multiple_metric_validation_model_fail(qc_extracted_metrics):
    """test MetricValidation for multiple metrics with failing conditions"""

    # GIVEN input attributes that does not meet the specified conditions
    metrics = copy.deepcopy(qc_extracted_metrics)
    metrics[4]["value"] = 2.0  # GC_DROPOUT set to 2.0 (failing condition)
    metrics[8]["value"] = 0.5  # PCT_TARGET_BASES_500X set to 50% (failing condition)

    # THEN check that the model filters the metrics according to its norm
    with pytest.raises(ValueError) as val_exc:
        MetricValidation(metrics=metrics)
    assert "2 validation errors for MetricValidation" in str(val_exc.value)
    assert metrics[4]["name"] in str(val_exc.value)
    assert metrics[8]["name"] in str(val_exc.value)


def test_metric_validation_model_norm_fail(qc_extracted_metrics):
    """test MetricValidation ValueError raising for an operator that it is not accepted"""

    # GIVEN a metric with an incorrect norm attribute
    metrics = copy.deepcopy(qc_extracted_metrics)
    metrics[4]["condition"]["norm"] = "lower"

    # THEN model raises an error due to a non accepted norm
    try:
        MetricValidation(metrics=metrics)
    except KeyError as key_exc:
        assert metrics[4]["condition"]["norm"] in str(key_exc)
