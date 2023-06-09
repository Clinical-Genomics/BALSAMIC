import copy
import json
import os.path
from pathlib import Path

import yaml

from BALSAMIC.assets.scripts.collect_qc_metrics import (
    get_multiqc_data_source,
    get_multiqc_metrics,
    collect_qc_metrics,
    get_qc_supported_capture_kit,
    get_requested_metrics,
    extract_number_variants,
    get_variant_metrics,
    get_metric_condition,
    get_relatedness_metrics,
)


def test_get_qc_supported_capture_kit(qc_requested_metrics):
    """test extraction of the capture kit name available for analysis"""

    # GIVEN a capture kit
    capture_kit = "panel_1_v1.0_hg19_design.bed"

    # GIVEN an expected output
    expected_output = "panel_1"

    # WHEN calling the function
    supported_capture_kit = get_qc_supported_capture_kit(
        capture_kit, qc_requested_metrics["targeted"]
    )

    # THEN check if the extracted bed file name corresponds to the expected one
    assert supported_capture_kit == expected_output


def test_get_requested_metrics_targeted(config_dict, qc_requested_metrics):
    """test retrieval of the requested targeted metrics"""

    # GIVEN a config_dict
    config = copy.deepcopy(config_dict)
    config["panel"]["capture_kit"] = "tests/panel/panel_2_v1.0_hg19_design.bed"

    # GIVEN the expected output
    expected_output = {
        "METRIC_1": {"condition": {"norm": "gt", "threshold": 1}},
        "METRIC_2": {"condition": {"norm": "gt", "threshold": 22}},
        "METRIC_4": {"condition": {"norm": "gt", "threshold": 4}},
    }

    # WHEN calling the function
    requested_metrics = get_requested_metrics(config, qc_requested_metrics)

    # THEN check if the requested targeted metrics are correctly retrieved
    assert requested_metrics.items() == expected_output.items()


def test_get_requested_metrics_wgs(config_dict, qc_requested_metrics):
    """test extraction of the requested WGS metrics"""

    # GIVEN a config_dict
    config = copy.deepcopy(config_dict)
    config["analysis"]["sequencing_type"] = "wgs"
    config["panel"] = None

    # GIVEN the expected output
    expected_output = {
        "METRIC_1": {"condition": {"norm": "gt", "threshold": 1}},
    }

    # WHEN calling the function
    requested_metrics = get_requested_metrics(config, qc_requested_metrics)

    # THEN check if the requested metrics are WGS specific
    assert requested_metrics.items() == expected_output.items()


def test_get_metric_condition(config_dict, qc_requested_metrics):
    """test metric condition extraction for WGS"""

    # GIVEN a config_dict
    config = copy.deepcopy(config_dict)
    seq_type = "wgs"
    config["analysis"]["sequencing_type"] = seq_type
    config["panel"] = None

    # GIVEN a specific sample & metric name
    sample = "ACC1"
    metric = "METRIC_1"

    # GIVEN an expected output
    expected_output = {"norm": "gt", "threshold": 1}

    # WHEN calling the function
    metric_condition = get_metric_condition(
        config, qc_requested_metrics[seq_type], sample, metric
    )

    # THEN check if the requested metrics has been correctly identified
    assert metric_condition.items() == expected_output.items()


def test_get_metric_condition_pct_wgs(config_dict):
    """test metric condition extraction for WGS"""

    # GIVEN a config_dict
    config = copy.deepcopy(config_dict)
    seq_type = "wgs"
    config["analysis"]["sequencing_type"] = seq_type
    config["panel"] = None

    # GIVEN a specific sample & metric name
    sample = "ACC1"
    metric = "PCT_60X"

    # GIVEN the requested metric
    req_metric = {
        "PCT_60X": {"condition": {"norm": "gt", "threshold": 1}},
    }

    # GIVEN an expected output
    expected_output = {"norm": "gt", "threshold": 1}

    # WHEN calling the function
    metric_condition = get_metric_condition(config, req_metric, sample, metric)

    # THEN check if the requested metrics has been correctly identified
    assert metric_condition.items() == expected_output.items()


def test_get_multiqc_data_source(multiqc_data_path):
    """test multiqc source extraction from multiqc_data.json analysis file"""

    # GIVEN input parameters and the multiqc data
    sample = "tumor.ACC1"
    source_name_hs_metrics = "multiqc_picard_HsMetrics"
    source_name_dup = "multiqc_picard_dups"

    with open(multiqc_data_path, "r") as f:
        multiqc_data = json.load(f)

    # GIVEN an expected output
    source_hs_metrics = "ACC1.dedup.realign.hsmetric"
    source_dup = "tumor.ACC1.dedup.metrics"

    # WHEN extracting the source of a specific sample and collection of metrics
    out_source_hs_metrics = get_multiqc_data_source(
        multiqc_data, sample, source_name_hs_metrics
    )
    out_source_dup = get_multiqc_data_source(multiqc_data, sample, source_name_dup)

    # THEN check if the extracted source names correspond to the expected ones
    assert source_hs_metrics == out_source_hs_metrics
    assert source_dup == out_source_dup

def test_get_multiqc_metrics(config_dict, multiqc_data_dict, qc_extracted_metrics):
    """test metrics retrieval from the multiqc_data.json file"""

    # GIVEN a config_dict
    config = copy.deepcopy(config_dict)
    config["panel"]["capture_kit"] = "tests/panel/lymphoma_6.1_hg19_design.bed"

    # WHEN calling the function
    metrics = get_multiqc_metrics(config, multiqc_data_dict)

    # THEN check if the metrics are correctly retrieved
    assert qc_extracted_metrics == metrics


def test_get_multiqc_metrics_filtering_umi(config_dict, multiqc_data_dict):
    """tests that UMI data is filtered out when extracting metrics"""

    # GIVEN a config_dict

    # WHEN calling the function
    metrics = get_multiqc_metrics(config_dict, multiqc_data_dict)

    # THEN check if the UMI samples are filtered out
    for metric in metrics:
        assert "umi" not in metric["input"]


def test_extract_number_variants():
    """tests number of variants formatting"""

    # GIVEN a raw input list of variant metrics
    counts = [
        "Number of samples: 2",
        "Number of SNPs:    111",
        "Number of INDELs:  14",
        "Number of MNPs:    0",
        "Number of sites:   125",
        "",
    ]

    # GIVEN an expected output after arranging the input list
    expected_variants_metrics = {
        "NUMBER_OF_SAMPLES": 2,
        "NUMBER_OF_SNPS": 111,
        "NUMBER_OF_INDELS": 14,
        "NUMBER_OF_MNPS": 0,
        "NUMBER_OF_SITES": 125,
    }

    # WHEN performing the extraction of variant metrics
    variant_metrics = extract_number_variants(counts)

    # THEN verify that the number of variants has been correctly retrieved
    assert expected_variants_metrics == variant_metrics


def test_get_variant_metrics(bcftools_counts_path):
    """tests variant metrics retrieval"""

    # GIVEN an SVDB bcftools counts path

    # GIVEN an expected MetricsModel dictionary
    expected_output_metris = {
        "header": None,
        "id": "case",
        "input": os.path.basename(bcftools_counts_path),
        "name": "NUMBER_OF_SITES",
        "step": "collect_custom_qc_metrics",
        "value": 125,
        "condition": {"norm": "lt", "threshold": 50000.0},
    }

    # WHEN extracting the number of variants
    output_metrics = get_variant_metrics(bcftools_counts_path)

    # THEN check that the output metrics has been correctly shaped
    assert expected_output_metris == output_metrics[0]


def test_collect_qc_metrics_targeted(
    tmp_path, config_path, multiqc_data_path, cli_runner
):
    """tests qc metrics yaml file generation for targeted analysis"""

    # GIVEN the output and multiqc metrics paths
    output_path = str(tmp_path / "sample_tumor_only_metrics_deliverables.yaml")

    # GIVEN a config path

    # WHEN invoking the python script
    result = cli_runner.invoke(
        collect_qc_metrics,
        [config_path, output_path, multiqc_data_path],
    )

    # THEN check if the YAML is correctly created and there are no errors
    assert result.exit_code == 0
    assert Path(output_path).exists()


def test_collect_qc_metrics_counts(
    tmp_path, config_path, multiqc_data_path, bcftools_counts_path, cli_runner
):
    """tests qc metrics yaml file generation for targeted analysis and providing a bcftools counts path"""

    # GIVEN the output, multiqc metrics and bcftools counts paths
    output_path = str(tmp_path / "sample_tumor_only_metrics_deliverables.yaml")

    # GIVEN a config path

    # WHEN invoking the python script
    result = cli_runner.invoke(
        collect_qc_metrics,
        [
            config_path,
            output_path,
            multiqc_data_path,
            bcftools_counts_path,  # multiple counts path regarding different variant callers
            bcftools_counts_path,
            bcftools_counts_path,
        ],
    )

    # THEN check if the YAML is correctly created and there are no errors
    assert result.exit_code == 0
    assert Path(output_path).exists()


def test_get_relatedness_metrics(multiqc_data_dict):
    """Tests relatedness metrics retrieval."""

    # GIVEN a multiqc_data_dict

    # GIVEN an expected MetricsModel dictionary
    expected_relatedness_metric = [
        {
            "header": None,
            "id": "id1",
            "input": "somalier.pairs.tsv",
            "name": "RELATEDNESS",
            "step": "multiqc_somalier",
            "value": 1,
            "condition": {"norm": "gt", "threshold": 0.80},
        }
    ]

    # WHEN extracting the relatedness metric
    relatedness_metric = get_relatedness_metrics(multiqc_data_dict)

    # THEN check that the relatedness metrics has been correctly shaped
    assert expected_relatedness_metric == relatedness_metric
