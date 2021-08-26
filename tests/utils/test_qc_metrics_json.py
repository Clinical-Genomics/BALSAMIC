import json

from BALSAMIC.utils.qc_metrics import get_qc_metrics_json


def test_get_qc_metrics_json():
    """test JSON object generation"""
    # GIVEN analysis_path and sequencing type
    analysis_path = "tests/test_data/qc_files/analysis"
    sequencing_type = "wgs"

    # WHEN invokes function
    qc_metrics = get_qc_metrics_json(analysis_path, sequencing_type)

    # THEN check if the obtained metrics are a valid JSON object
    try:
        json.loads(qc_metrics)
        assert True
    except TypeError:
        assert False
