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
    capture_kit_resolve_type,
    extract_number_variants,
    get_variant_metrics,
)


def test_capture_kit_resolve_type():
    """test capture_kit type"""

    # GIVEN an expected output
    capture_kit = "panel.bed"

    # THEN check if the extracted capture kit is correctly formatted
    assert capture_kit_resolve_type("None") is None
    assert capture_kit_resolve_type(capture_kit) == capture_kit


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


def test_get_requested_metrics_targeted(qc_requested_metrics):
    """test retrieval of the requested targeted metrics"""

    # GIVEN a sequencing type and a capture kit
    seq_type = "targeted"
    capture_kit = "panel_2_v1.0_hg19_design.bed"

    # GIVEN the expected output
    expected_output = {
        "METRIC_1": {"condition": {"norm": "gt", "threshold": 1}},
        "METRIC_2": {"condition": {"norm": "gt", "threshold": 22}},
        "METRIC_4": {"condition": {"norm": "gt", "threshold": 4}},
    }

    # WHEN calling the function
    requested_metrics = get_requested_metrics(
        qc_requested_metrics, seq_type, capture_kit
    )

    # THEN check if the requested targeted metrics are correctly retrieved
    assert requested_metrics.items() == expected_output.items()


def test_get_requested_metrics_wgs(qc_requested_metrics):
    """test extraction of the requested WGS metrics"""

    # GIVEN a sequencing type and a capture kit
    seq_type = "wgs"
    capture_kit = None

    # GIVEN the expected output
    expected_output = {
        "METRIC_1": {"condition": {"norm": "gt", "threshold": 1}},
    }

    # WHEN calling the function
    requested_metrics = get_requested_metrics(
        qc_requested_metrics, seq_type, capture_kit
    )

    # THEN check if the requested metrics are WGS specific
    assert requested_metrics.items() == expected_output.items()


def test_get_multiqc_data_source(multiqc_data_path):
    """test multiqc source extraction from multiqc_data.json analysis file"""

    # GIVEN input parameters and the multiqc data
    sample = "concatenated_tumor_XXXXXX_R"
    source_name_hs_metrics = "multiqc_picard_HsMetrics"
    source_name_dup = "multiqc_picard_dups"

    with open(multiqc_data_path, "r") as f:
        multiqc_data = json.load(f)

    # GIVEN an expected output
    source_hs_metrics = "concatenated_tumor_XXXXXX_R.sorted.mrkdup.hsmetric"
    source_dup = "concatenated_tumor_XXXXXX_R.sorted.mrkdup.txt"

    # WHEN extracting the source of a specific sample and collection of metrics
    out_source_hs_metrics = get_multiqc_data_source(
        multiqc_data, sample, source_name_hs_metrics
    )
    out_source_dup = get_multiqc_data_source(multiqc_data, sample, source_name_dup)

    # THEN check if the extracted source names correspond to the expected ones
    assert source_hs_metrics == out_source_hs_metrics
    assert source_dup == out_source_dup


def test_get_multiqc_metrics(multiqc_data_path, qc_extracted_metrics):
    """test metrics retrieval from the multiqc_data.json file"""

    # GIVEN a sequencing type and a capture kit
    seq_type = "targeted"
    capture_kit = "lymphoma_6.1_hg19_design.bed"

    # WHEN calling the function
    metrics = get_multiqc_metrics(
        multiqc_data_path,
        seq_type,
        capture_kit,
    )

    # THEN check if the metrics are correctly retrieved
    assert qc_extracted_metrics == metrics


def test_get_multiqc_metrics_filtering_umi(multiqc_data_path):
    """tests that UMI data is filtered out when extracting metrics"""

    # GIVEN a sequencing type and a capture kit
    seq_type = "targeted"
    capture_kit = None

    # WHEN calling the function
    metrics = get_multiqc_metrics(
        multiqc_data_path,
        seq_type,
        capture_kit,
    )

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


def test_collect_qc_metrics_targeted(tmp_path, multiqc_data_path, cli_runner):
    """tests qc metrics yaml file generation for targeted analysis"""

    # GIVEN the output and multiqc metrics paths
    output_path = tmp_path / "sample_tumor_only_metrics_deliverables.yaml"

    # GIVEN a sequencing type and a capture kit
    seq_type = "targeted"
    capture_kit = "lymphoma_6.1_hg19_design.bed"

    # WHEN invoking the python script
    result = cli_runner.invoke(
        collect_qc_metrics,
        [str(output_path), multiqc_data_path, seq_type, capture_kit],
    )

    # THEN check if the YAML is correctly created and there are no errors
    assert result.exit_code == 0
    assert Path(output_path).exists()


def test_collect_qc_metrics_wgs(tmp_path, multiqc_data_path, cli_runner):
    """tests qc metrics yaml file generation for wgs analysis"""

    # GIVEN the output and multiqc metrics paths
    output_path = tmp_path / "sample_tumor_only_wgs_metrics_deliverables.yaml"

    # GIVEN a sequencing type and a capture kit
    seq_type = "wgs"
    capture_kit = "None"

    # WHEN invoking the python script
    result = cli_runner.invoke(
        collect_qc_metrics,
        [str(output_path), multiqc_data_path, seq_type, capture_kit],
    )

    # THEN check if the YAML is correctly created and there are no errors
    assert result.exit_code == 0
    assert Path(output_path).exists()


def test_collect_qc_metrics_counts(
    tmp_path, multiqc_data_path, bcftools_counts_path, cli_runner
):
    """tests qc metrics yaml file generation for targeted analysis and providing a bcftools counts path"""

    # GIVEN the output, multiqc metrics and bcftools counts paths
    output_path = tmp_path / "sample_tumor_only_metrics_deliverables.yaml"

    # GIVEN a sequencing type and a capture kit
    seq_type = "targeted"
    capture_kit = "gmsmyeloid_5.2_hg19_design.bed"

    # WHEN invoking the python script
    result = cli_runner.invoke(
        collect_qc_metrics,
        [
            str(output_path),
            multiqc_data_path,
            bcftools_counts_path,  # multiple counts path regarding different variant callers
            bcftools_counts_path,
            bcftools_counts_path,
            seq_type,
            capture_kit,
        ],
    )

    # THEN check if the YAML is correctly created and there are no errors
    assert result.exit_code == 0
    assert Path(output_path).exists()
