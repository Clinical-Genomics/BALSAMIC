import json
import os

from BALSAMIC.utils.qc_metrics import (
    get_qc_metrics_json,
    read_metrics,
    update_metrics_dict,
    get_qc_metrics_dict,
    get_multiqc_data_source,
    extract_metrics_for_delivery,
)


def test_read_metrics(analysis_path):
    """test metric extraction from a specific QC file"""

    # GIVEN a QC file name
    file_name = "multiqc_picard_dups.json"

    # Given an expected output
    expected_output = {
        "concatenated_tumor_XXXXXX_R": {
            "LIBRARY": "Unknown Library",
            "UNPAIRED_READS_EXAMINED": 11860.0,
            "READ_PAIRS_EXAMINED": 20440841.0,
            "SECONDARY_OR_SUPPLEMENTARY_RDS": 4333388.0,
            "UNMAPPED_READS": 19824.0,
            "UNPAIRED_READ_DUPLICATES": 10178.0,
            "READ_PAIR_DUPLICATES": 14680829.0,
            "READ_PAIR_OPTICAL_DUPLICATES": 0.0,
            "PERCENT_DUPLICATION": 0.718251,
            "ESTIMATED_LIBRARY_SIZE": 5951948.0,
        }
    }

    # WHEN calling the function
    raw_metrics = read_metrics(analysis_path, file_name)

    # THEN check if the extracted metrics correspond to the expected ones
    assert raw_metrics.items() == expected_output.items()


def test_update_metrics_dict(qc_extracted_metrics):
    """test adding metrics to a nested dictionary"""

    # GIVEN input parameters
    sample_id = "sample_"
    metric = ["MEAN_INSERT_SIZE", {"condition": {"norm": "lt", "threshold": 1.0}}]
    value = 0.5

    # WHEN adding a metric to an empty dictionary
    metric[0] = "MEAN_INSERT_SIZE_1"
    m_dict = update_metrics_dict(sample_id + "1", metric, value, {})

    # WHEN appending a metric to an already created dictionary
    metric[0] = "MEAN_INSERT_SIZE_2"
    m_dict = update_metrics_dict(sample_id + "1", metric, value, m_dict)

    # WHEN appending a metric from another sample to a dictionary
    metric[0] = "MEAN_INSERT_SIZE_1"
    m_dict = update_metrics_dict(sample_id + "2", metric, value, m_dict)

    # THEN check if the dictionary is updated correctly
    assert m_dict.items() == qc_extracted_metrics["metrics"].items()


def test_get_qc_metrics_dict(analysis_path, qc_metrics):
    """test QC metric extraction and its format"""

    # GIVEN a sequencing type
    seq_type = "targeted"

    # GIVEN an expected output
    expected_output = {
        "concatenated_tumor": [
            {
                "name": "MEAN_INSERT_SIZE",
                "value": 74.182602,
                "condition": None,
                "meets_condition": None,
            },
            {
                "name": "MEAN_TARGET_COVERAGE",
                "value": 832.13854,
                "condition": {"norm": "gt", "threshold": 500.0},
                "meets_condition": None,
            },
        ]
    }

    # WHEN calling the function
    metrics_dict = get_qc_metrics_dict(analysis_path, qc_metrics["qc"][seq_type])

    # THEN check if the extracted metrics and its structure meets the expected one
    assert metrics_dict.items() == expected_output.items()


def test_get_qc_metrics_json(analysis_path):
    """test JSON object generation"""

    # GIVEN a sequencing type
    seq_type = "targeted"

    # WHEN calling the function
    qc_metrics = get_qc_metrics_json(analysis_path, seq_type)

    # THEN check if the obtained metrics are a valid JSON object
    try:
        json.loads(qc_metrics)
        assert True
    except TypeError:
        assert False


def test_get_multiqc_data_source(analysis_path):
    """test multiQC source extraction from multiqc_data.json analysis file"""

    # GIVEN input parameters
    sample = "concatenated_normal_XXXXXX_R"
    source_name_hs_metrics = "multiqc_picard_HsMetrics"
    source_name_dup = "multiqc_picard_dups"

    with open(
        os.path.join(analysis_path, "qc", "multiqc_data", "multiqc_data.json"), "r"
    ) as f:
        raw_data = json.load(f)

    # GIVEN an expected output
    source_hs_metrics = "concatenated_normal_XXXXXX_R.umi.collect_hsmetric"
    source_dup = "concatenated_normal_XXXXXX_R.sorted.mrkdup.txt"

    # WHEN extracting the source of a specific sample and collection of metrics
    out_source_hs_metrics = get_multiqc_data_source(
        raw_data, sample, source_name_hs_metrics
    )
    data_source_dup = get_multiqc_data_source(raw_data, sample, source_name_dup)

    # THEN check if the extracted source names correspond to the expected ones
    assert source_hs_metrics == out_source_hs_metrics
    assert source_dup == data_source_dup


def test_extract_metrics_for_delivery(analysis_path):
    """test output metrics retrieving"""

    # GIVEN a sequencing type
    seq_type = "targeted"

    # GIVEN an expected output
    n_metrics = 22

    first_metric = {
        "header": None,
        "id": "normal",
        "input": "concatenated_normal_XXXXXX_R.umi.collect_hsmetric",
        "name": "PCT_OFF_BAIT",
        "step": "multiqc_picard_HsMetrics",
        "value": 0.401673,
    }

    last_metric = {
        "header": None,
        "id": "tumor",
        "input": "concatenated_tumor_XXXXXX_R.sorted.mrkdup.txt",
        "name": "PERCENT_DUPLICATION",
        "step": "multiqc_picard_dups",
        "value": 0.931044,
    }

    # WHEN calling the function
    metrics = extract_metrics_for_delivery(analysis_path, seq_type)

    # THEN check if the metrics are correctly retrieved
    assert len(metrics) == n_metrics
    assert first_metric in metrics and last_metric in metrics
