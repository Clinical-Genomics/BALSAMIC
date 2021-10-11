import json

from BALSAMIC.utils.qc_metrics import (
    get_qc_metrics_json,
    read_metrics,
    update_metrics_dict,
    get_qc_metrics_dict,
    get_qc_filtered_metrics_json,
    merge_dicts,
    get_qc_available_panel_beds,
)


def test_get_qc_available_panel_beds():
    """test extraction of the panel beds available for QC validation"""

    # GIVEN quality control metrics
    metrics = {
        "common": {
            "metrics_1.json": {"METRIC_1": 0.1, "METRIC_2": 0.2},
        },
        "panel_1.bed": {"metrics_1.json": {"METRIC_4": 0.4}},
        "panel_2.bed": {"metrics_3.json": {"METRIC_5": 0.5}},
    }

    # GIVEN an expected output
    expected_output = ["panel_1.bed", "panel_2.bed"]

    # WHEN calling the function
    available_panel_beds = get_qc_available_panel_beds(metrics)

    # THEN check if the extracted bed file names correspond to the expected ones
    assert available_panel_beds == expected_output


def test_read_metrics(analysis_path):
    """test metric extraction from a specific QC file"""

    # GIVEN a QC file name
    file_name = "multiqc_picard_dups.json"

    # GIVEN an expected output
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


def test_merge_dicts():
    """test dictionary merging"""

    # GIVEN input dictionary
    nested_dict = {
        "common": {
            "metrics_1.json": {"METRIC_1": 0.1, "METRIC_2": 0.2},
            "metrics_2.json": {"METRIC_3": 0.3},
        },
        "panel_1.bed": {"metrics_2.json": {"METRIC_4": 0.4}},
        "panel_2.bed": {  # Tests that the common metric condition is overwritten
            "metrics_1.json": {"METRIC_1": 0.5}
        },
    }

    # GIVEN an expected output
    expected_output = {
        "metrics_1.json": {"METRIC_1": 0.5, "METRIC_2": 0.2},
        "metrics_2.json": {"METRIC_3": 0.3, "METRIC_4": 0.4},
    }

    # WHEN calling the function
    merged_dict = merge_dicts(
        nested_dict["common"], nested_dict["panel_1.bed"], nested_dict["panel_2.bed"]
    )

    # THEN check if the extracted output meets the merged dictionary
    assert merged_dict.items() == expected_output.items()


def test_get_qc_metrics_json_targeted(analysis_path):
    """test JSON object generation"""

    # GIVEN a sequencing type
    seq_type = "targeted"
    capture_kit = "gicfdna_3.1_hg19_design.bed"

    # WHEN calling the function
    qc_metrics = get_qc_metrics_json(
        analysis_path, seq_type, capture_kit
    )

    # THEN check if the obtained metrics are a valid JSON object
    try:
        assert sorted(
            json.loads(get_qc_filtered_metrics_json(qc_metrics, "failed"))[
                "concatenated_tumor"
            ].keys()
        ) == sorted(["FOLD_80_BASE_PENALTY", "PCT_OFF_BAIT"])
    except TypeError:
        assert False


def test_get_qc_metrics_json_WGS(analysis_path):
    """test JSON object generation for a custom bed"""

    # GIVEN a sequencing type
    seq_type = "wgs"
    capture_kit = None

    # WHEN calling the function
    qc_metrics = get_qc_metrics_json(analysis_path, seq_type, capture_kit)

    # THEN check if the obtained metrics are WGS specific
    try:
        assert sorted(
            json.loads(get_qc_filtered_metrics_json(qc_metrics, "passed"))[
                "concatenated_tumor"
            ].keys()
        ) == sorted(["MEAN_INSERT_SIZE", "PERCENT_DUPLICATION"])
    except TypeError:
        assert False


def test_get_qc_filtered_metrics():
    """test filtered metric extraction from an already generated JSON object"""

    # GIVEN a JSON object
    qc_metrics_json = json.dumps(
        {
            "sample_1": {"failed": {}, "passed": {"METRIC_1": 0.5}},
            "sample_2": {
                "failed": {"METRIC_1": 0.5},
                "passed": {"METRIC_2": 0.5, "METRIC_3": 0.5},
            },
        }
    )

    # GIVEN the expected output
    expected_output = {
        "sample_1": {"METRIC_1": 0.5},
        "sample_2": {"METRIC_2": 0.5, "METRIC_3": 0.5},
    }

    # WHEN calling the function
    qc_passed_metrics = get_qc_filtered_metrics_json(qc_metrics_json, "passed")

    # THEN check if the obtained metrics are the ones that passed the QC validation
    assert json.loads(qc_passed_metrics).items() == expected_output.items()


def test_get_qc_filtered_metrics_empty():
    """test empty return when extracting filtered metrics"""

    # GIVEN a JSON object
    qc_metrics_json = json.dumps(
        {
            "sample_1": {"failed": {}, "passed": {"METRIC_1": 0.5, "METRIC_2": 0.5}},
            "sample_2": {"failed": {}, "passed": {"METRIC_2": 0.5}},
        }
    )

    # GIVEN the expected output
    expected_output = {}

    # WHEN calling the function
    qc_failed_metrics = get_qc_filtered_metrics_json(qc_metrics_json, "failed")

    # THEN check if the are no output metrics
    assert not json.loads(qc_failed_metrics)
    assert json.loads(qc_failed_metrics).items() == expected_output.items()
