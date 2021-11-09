import json
import os

from pydantic import ValidationError

from BALSAMIC.utils.qc_metrics import (
    get_qc_metrics_json,
    read_metrics,
    update_metrics_dict,
    get_qc_metrics_dict,
    get_qc_available_panel_beds,
    merge_dicts,
    get_multiqc_data_source,
    extract_metrics_for_delivery,
)


def test_get_qc_available_panel_beds(qc_raw_targeted_metrics):
    """test extraction of the panel beds available for QC validation"""

    # GIVEN an expected output
    expected_output = ["panel_1.bed", "panel_2.bed"]

    # WHEN calling the function
    available_panel_beds = get_qc_available_panel_beds(qc_raw_targeted_metrics)

    # THEN check if the extracted bed file names correspond to the expected ones
    assert available_panel_beds == expected_output


def test_merge_dicts(qc_raw_targeted_metrics):
    """test dictionary merging and requirements overwriting by panel BED specific conditions"""

    # GIVEN an expected output
    expected_output = {
        "metrics_1.json": {"METRIC_1": 0.5, "METRIC_2": 0.2, "METRIC_4": 0.4},
        "metrics_2.json": {"METRIC_3": 0.3},
    }

    # WHEN calling the function
    merged_dict = merge_dicts(
        qc_raw_targeted_metrics["default"],
        qc_raw_targeted_metrics["panel_2.bed"],
    )

    # THEN check if the extracted output meets the merged dictionary
    assert merged_dict.items() == expected_output.items()


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
    """test QC metric extraction and its structure"""

    # GIVEN a sequencing type
    seq_type = "targeted"

    # GIVEN an expected output
    expected_output = {
        "concatenated_tumor": [
            {
                "name": "MEAN_INSERT_SIZE",
                "norm": None,
                "threshold": None,
                "value": 74.182602,
            },
            {
                "name": "MEDIAN_TARGET_COVERAGE",
                "norm": "gt",
                "threshold": 500.0,
                "value": 461.0,
            },
        ]
    }

    # WHEN calling the function
    metrics_dict = get_qc_metrics_dict(analysis_path, qc_metrics["qc"][seq_type])

    # THEN check if the extracted metrics and its structure meets the expected one
    assert metrics_dict.items() == expected_output.items()


def test_get_qc_metrics_json_wgs(analysis_path):
    """test JSON object generation for a WGS run"""

    # GIVEN a sequencing type
    seq_type = "wgs"
    capture_kit = None

    # GIVEN retrieved WGS metrics
    output_metrics = {"concatenated_tumor": {"FOLD_80_BASE_PENALTY": 1.238604}}

    # WHEN calling the function
    qc_metrics = get_qc_metrics_json(analysis_path, seq_type, capture_kit)

    # THEN check if the obtained metrics are WGS specific
    assert qc_metrics.items() == output_metrics.items()


def test_get_qc_metrics_json_targeted(analysis_path):
    """test JSON object generation for a custom bed file"""

    # GIVEN a sequencing type
    seq_type = "targeted"
    capture_kit = "lymphoma_6.1_hg19_design.bed"

    # THEN check if the obtained metrics are following the panel bed specific requirements
    try:
        get_qc_metrics_json(analysis_path, seq_type, capture_kit)
    except ValidationError as val_err:
        assert (
            "2 validation errors for QCValidationModel" in str(val_err)
            and "MEDIAN_TARGET_COVERAGE" in str(val_err)
            and "FOLD_80_BASE_PENALTY" in str(val_err)
        )


def test_get_multiqc_data_source(analysis_path):
    """test multiQC source extraction from multiqc_data.json analysis file"""

    # GIVEN input parameters
    sample = "concatenated_tumor_XXXXXX_R"
    source_name_hs_metrics = "multiqc_picard_HsMetrics"
    source_name_dup = "multiqc_picard_dups"

    with open(
        os.path.join(analysis_path, "qc", "multiqc_data", "multiqc_data.json"), "r"
    ) as f:
        raw_data = json.load(f)

    # GIVEN an expected output
    source_hs_metrics = "concatenated_tumor_XXXXXX_R.sorted.mrkdup.hsmetric"
    source_dup = "concatenated_tumor_XXXXXX_R.sorted.mrkdup.txt"

    # WHEN extracting the source of a specific sample and collection of metrics
    out_source_hs_metrics = get_multiqc_data_source(
        raw_data, sample, source_name_hs_metrics
    )
    out_source_dup = get_multiqc_data_source(raw_data, sample, source_name_dup)

    # THEN check if the extracted source names correspond to the expected ones
    assert source_hs_metrics == out_source_hs_metrics
    assert source_dup == out_source_dup


def test_extract_metrics_for_delivery(analysis_path):
    """test output metrics retrieving"""

    # GIVEN a sequencing type
    seq_type = "targeted"

    # GIVEN an expected output
    n_metrics = 6  # Number of expected metric

    hs_metric = {
        "header": None,
        "id": "tumor",
        "input": "concatenated_tumor_XXXXXX_R.sorted.mrkdup.hsmetric",
        "name": "PCT_OFF_BAIT",
        "step": "multiqc_picard_HsMetrics",
        "value": 0.364546,
    }

    ins_size_metric = {
        "header": None,
        "id": "tumor",
        "input": "concatenated_tumor_XXXXXX_R.sorted.insertsizemetric",
        "name": "MEAN_INSERT_SIZE",
        "step": "multiqc_picard_insertSize",
        "value": 201.813054,
    }

    dups_metric = {
        "header": None,
        "id": "tumor",
        "input": "concatenated_tumor_XXXXXX_R.sorted.mrkdup.txt",
        "name": "PERCENT_DUPLICATION",
        "step": "multiqc_picard_dups",
        "value": 0.391429,
    }

    # WHEN calling the function
    metrics = extract_metrics_for_delivery(analysis_path, seq_type)

    # THEN check if the metrics are correctly retrieved
    assert len(metrics) == n_metrics
    assert (
        hs_metric in metrics and ins_size_metric in metrics and dups_metric in metrics
    )


def test_extract_metrics_for_delivery_filtering_umi(analysis_path):
    """test umi discarding when extracting metrics"""

    # GIVEN a sequencing type
    seq_type = "targeted"

    # WHEN calling the function
    metrics = extract_metrics_for_delivery(analysis_path, seq_type)

    # THEN check if the umi samples are filtered out
    for metric in metrics:
        assert "umi" not in metric["input"]
