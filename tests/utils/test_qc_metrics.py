from BALSAMIC.utils.qc_metrics import validate_qc_metrics


def test_validate_qc_metrics(qc_extracted_metrics):
    """test QC metric validation"""

    # WHEN calling the function
    validated_metrics_pass = validate_qc_metrics(qc_extracted_metrics)

    # THEN check if the obtained metrics are correctly parsed and validated
    assert validated_metrics_pass == qc_extracted_metrics
