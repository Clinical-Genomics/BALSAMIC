import json
from pathlib import Path
from unittest import mock

import graphviz

from BALSAMIC.constants.constants import FileType
from tests.conftest import MOCKED_OS_ENVIRON


def test_tumor_normal_config(
    invoke_cli,
    case_id_tumor_normal: str,
    tumor_sample_name: str,
    normal_sample_name: str,
    analysis_dir: str,
    fastq_dir_tumor_normal_parameterize: str,
    balsamic_cache: str,
    panel_bed_file: str,
    sentieon_license: str,
    sentieon_install_dir: str,
):
    """Test tumor normal balsamic config case command."""

    # GIVEN a case ID, fastq files, and an analysis dir
    # WHEN creating a case analysis
    with mock.patch.dict(
        MOCKED_OS_ENVIRON,
        {
            "SENTIEON_LICENSE": sentieon_license,
            "SENTIEON_INSTALL_DIR": sentieon_install_dir,
        },
    ):
        result = invoke_cli(
            [
                "config",
                "case",
                "--case-id",
                case_id_tumor_normal,
                "--gender",
                "male",
                "--analysis-dir",
                analysis_dir,
                "--fastq-path",
                fastq_dir_tumor_normal_parameterize,
                "-p",
                panel_bed_file,
                "--balsamic-cache",
                balsamic_cache,
                "--tumor-sample-name",
                tumor_sample_name,
                "--normal-sample-name",
                normal_sample_name,
            ],
        )

    # THEN a config should be created and exist
    assert result.exit_code == 0
    assert Path(
        analysis_dir, case_id_tumor_normal, f"{case_id_tumor_normal}.{FileType.JSON}"
    ).exists()


def test_tumor_only_config(
    invoke_cli,
    case_id_tumor_only: str,
    tumor_sample_name: str,
    analysis_dir: str,
    fastq_dir_tumor_only: str,
    balsamic_cache: str,
    panel_bed_file: str,
    sentieon_license: str,
    sentieon_install_dir: str,
):
    """Test tumor only balsamic config case command."""

    # GIVEN a case ID, fastq files, and an analysis dir

    # WHEN creating a case analysis
    with mock.patch.dict(
        MOCKED_OS_ENVIRON,
        {
            "SENTIEON_LICENSE": sentieon_license,
            "SENTIEON_INSTALL_DIR": sentieon_install_dir,
        },
    ):
        result = invoke_cli(
            [
                "config",
                "case",
                "--case-id",
                case_id_tumor_only,
                "--analysis-dir",
                analysis_dir,
                "--fastq-path",
                fastq_dir_tumor_only,
                "-p",
                panel_bed_file,
                "--balsamic-cache",
                balsamic_cache,
                "--tumor-sample-name",
                tumor_sample_name,
            ],
        )

    # THEN a config should be created and exist
    assert result.exit_code == 0
    assert Path(
        analysis_dir, case_id_tumor_only, f"{case_id_tumor_only}.{FileType.JSON}"
    ).exists()


def test_run_without_permissions(
    invoke_cli,
    case_id_tumor_only: str,
    tumor_sample_name: str,
    fastq_dir_tumor_only: str,
    no_write_perm_path: str,
    panel_bed_file: str,
    balsamic_cache: str,
):
    """Test balsamic config case with no write permissions to the analysis directory."""

    # GIVEN CLI arguments including an analysis_dir without write permissions

    # WHEN invoking the config case command
    result = invoke_cli(
        [
            "config",
            "case",
            "--case-id",
            case_id_tumor_only,
            "--analysis-dir",
            no_write_perm_path,
            "--fastq-path",
            fastq_dir_tumor_only,
            "-p",
            panel_bed_file,
            "--balsamic-cache",
            balsamic_cache,
            "--tumor-sample-name",
            tumor_sample_name,
        ],
    )

    # THEN program exits before completion
    assert result.exit_code == 1


def test_tumor_only_umi_config_background_file(
    invoke_cli,
    case_id_tumor_only: str,
    tumor_sample_name: str,
    analysis_dir: str,
    fastq_dir_tumor_only: str,
    balsamic_cache: str,
    panel_bed_file: str,
    background_variant_file: str,
):
    """Test balsamic umi config case providing a background variants file."""

    # GIVEN CLI arguments including a background variant file

    # WHEN invoking the config case command
    result = invoke_cli(
        [
            "config",
            "case",
            "--case-id",
            case_id_tumor_only,
            "--analysis-workflow",
            "balsamic-umi",
            "--analysis-dir",
            analysis_dir,
            "--fastq-path",
            fastq_dir_tumor_only,
            "-p",
            panel_bed_file,
            "--background-variants",
            background_variant_file,
            "--balsamic-cache",
            balsamic_cache,
            "--tumor-sample-name",
            tumor_sample_name,
        ],
    )

    # THEN program exits and checks for filepath
    assert result.exit_code == 0
    assert Path(background_variant_file).exists()


def test_pon_cnn_file(
    invoke_cli,
    tumor_sample_name: str,
    analysis_dir: str,
    balsamic_cache: str,
    panel_bed_file: str,
    pon_cnn_path: str,
    fastq_dir_tumor_only: str,
    case_id_tumor_only: str,
):
    """Test balsamic config case with a PON reference."""

    # GIVEN CLI arguments including optional pon reference ".cnn" file

    # WHEN invoking the config case command
    result = invoke_cli(
        [
            "config",
            "case",
            "--case-id",
            case_id_tumor_only,
            "--analysis-dir",
            analysis_dir,
            "--fastq-path",
            fastq_dir_tumor_only,
            "-p",
            panel_bed_file,
            "--pon-cnn",
            pon_cnn_path,
            "--balsamic-cache",
            balsamic_cache,
            "--tumor-sample-name",
            tumor_sample_name,
        ],
    )

    # THEN program exits and checks for filepath
    assert result.exit_code == 0
    assert Path(pon_cnn_path).exists()


def test_dag_graph_success_tumor_only(tumor_only_config: str):
    """Test DAG graph building success."""
    # WHEN creating config using standard CLI input and setting Sentieon env vars

    # THEN DAG graph should be created successfully
    assert Path(json.load(open(tumor_only_config))["analysis"]["dag"]).exists()


def test_dag_graph_success_tumor_only_pon(tumor_only_pon_config: str):
    """Test DAG graph building success."""

    # WHEN creating config using standard CLI input and setting Sentieon env vars

    # THEN DAG graph should be created successfully
    assert Path(json.load(open(tumor_only_pon_config))["analysis"]["dag"]).exists()


def test_dag_graph_success_tumor_normal(tumor_normal_config: str):
    """Test DAG graph building success."""
    # WHEN creating config using standard CLI input and setting Sentieon env vars

    # THEN DAG graph should be created successfully
    assert Path(json.load(open(tumor_normal_config))["analysis"]["dag"]).exists()


def test_dag_graph_success_tumor_only_umi(tumor_only_umi_config: str):
    """Test DAG graph building success."""
    # WHEN creating config using standard CLI input and setting Sentieon env vars

    # THEN DAG graph should be created successfully
    assert Path(json.load(open(tumor_only_umi_config))["analysis"]["dag"]).exists()


def test_dag_graph_success_tumor_normal_umi(tumor_normal_umi_config: str):
    """Test DAG graph building success."""
    # WHEN creating config using standard CLI input and setting Sentieon env vars

    # THEN DAG graph should be created successfully
    assert Path(json.load(open(tumor_normal_umi_config))["analysis"]["dag"]).exists()


def test_dag_graph_success_tumor_only_wgs(
    tumor_only_wgs_config: str,
):
    """Test DAG graph building success."""
    # WHEN creating config using standard CLI input and setting Sentieon env vars

    # THEN DAG graph should be created successfully
    assert Path(json.load(open(tumor_only_wgs_config))["analysis"]["dag"]).exists()


def test_dag_graph_success_tumor_normal_wgs(tumor_normal_wgs_config: str):
    """Test DAG graph building success."""
    # WHEN creating config using standard CLI input and setting Sentieon env vars

    # THEN DAG graph should be created successfully
    assert Path(json.load(open(tumor_normal_wgs_config))["analysis"]["dag"]).exists()


def test_config_graph_failed(
    invoke_cli,
    case_id_tumor_only: str,
    tumor_sample_name: str,
    analysis_dir: str,
    fastq_dir_tumor_only: str,
    balsamic_cache: str,
    panel_bed_file: str,
):
    """Test DAG graph building failure."""

    # GIVEN an analysis config

    # GIVEN an empty graphviz instance
    with mock.patch.object(graphviz, "Source") as mocked:
        mocked.return_value = None
        case_result = invoke_cli(
            [
                "config",
                "case",
                "--case-id",
                case_id_tumor_only,
                "--analysis-dir",
                analysis_dir,
                "--fastq-path",
                fastq_dir_tumor_only,
                "-p",
                panel_bed_file,
                "--balsamic-cache",
                balsamic_cache,
                "--tumor-sample-name",
                tumor_sample_name,
            ],
        )

    # THEN the graph should not have been built
    assert case_result.exit_code == 1


def test_missing_required_gens_arguments(
    invoke_cli,
    tumor_sample_name: str,
    analysis_dir: str,
    balsamic_cache: str,
    fastq_dir_tumor_only: str,
    case_id_tumor_only: str,
    gens_cov_pon_file: str,
    gens_min_5_af_gnomad_file: str,
):
    """Test balsamic config case with 2 out of 3 required GENS arguments."""

    # GIVEN CLI arguments including optional GENS input-files

    # WHEN invoking the config case command
    result = invoke_cli(
        [
            "config",
            "case",
            "--case-id",
            case_id_tumor_only,
            "--analysis-dir",
            analysis_dir,
            "--fastq-path",
            fastq_dir_tumor_only,
            "--balsamic-cache",
            balsamic_cache,
            "--tumor-sample-name",
            tumor_sample_name,
            "--gens-coverage-pon",
            gens_cov_pon_file,
            "--gnomad-min-af5",
            gens_min_5_af_gnomad_file,
        ],
    )
    # THEN the CLI should exit code 2 and display an informative error message
    assert result.exit_code == 2
    assert (
        "All three arguments (genome_interval gens_coverage_pon, gnomad_min_af5) are required for GENS."
        in result.output
    )


def test_config_with_gens_arguments(
    invoke_cli,
    tumor_sample_name: str,
    analysis_dir: str,
    balsamic_cache: str,
    fastq_dir_tumor_only: str,
    case_id_tumor_only: str,
    gens_cov_pon_file: str,
    gens_min_5_af_gnomad_file: str,
    gens_hg19_interval_list: str,
):
    """Test balsamic config case with GENS arguments."""

    # GIVEN CLI arguments including optional GENS input-files

    # WHEN invoking the config case command
    result = invoke_cli(
        [
            "config",
            "case",
            "--case-id",
            case_id_tumor_only,
            "--analysis-dir",
            analysis_dir,
            "--fastq-path",
            fastq_dir_tumor_only,
            "--balsamic-cache",
            balsamic_cache,
            "--tumor-sample-name",
            tumor_sample_name,
            "--gens-coverage-pon",
            gens_cov_pon_file,
            "--gnomad-min-af5",
            gens_min_5_af_gnomad_file,
            "--genome-interval",
            gens_hg19_interval_list,
        ],
    )
    # THEN a config should be created and exist
    assert result.exit_code == 0
    assert Path(
        analysis_dir, case_id_tumor_only, f"{case_id_tumor_only}.{FileType.JSON}"
    ).exists()


def test_config_with_gens_arguments_for_tga(
    invoke_cli,
    tumor_sample_name: str,
    analysis_dir: str,
    balsamic_cache: str,
    fastq_dir_tumor_only: str,
    case_id_tumor_only: str,
    gens_cov_pon_file: str,
    gens_min_5_af_gnomad_file: str,
    gens_hg19_interval_list: str,
    panel_bed_file: str,
):
    """Test balsamic config case with GENS arguments for TGA."""

    # GIVEN CLI arguments including optional GENS input-files

    # WHEN invoking the config case command
    result = invoke_cli(
        [
            "config",
            "case",
            "--case-id",
            case_id_tumor_only,
            "--analysis-dir",
            analysis_dir,
            "--fastq-path",
            fastq_dir_tumor_only,
            "--balsamic-cache",
            balsamic_cache,
            "--tumor-sample-name",
            tumor_sample_name,
            "--gens-coverage-pon",
            gens_cov_pon_file,
            "--gnomad-min-af5",
            gens_min_5_af_gnomad_file,
            "--genome-interval",
            gens_hg19_interval_list,
            "-p",
            panel_bed_file,
        ],
    )
    # THEN config should fail with error message
    assert result.exit_code == 2
    assert (
        "GENS is currently not compatible with TGA analysis, only WGS." in result.output
    )


def test_config_wgs_with_exome(
    invoke_cli,
    tumor_sample_name: str,
    analysis_dir: str,
    balsamic_cache: str,
    fastq_dir_tumor_only: str,
    case_id_tumor_only: str,
):
    """Test balsamic config case with --exome argument for WGS."""

    # GIVEN CLI arguments including optional GENS input-files

    # WHEN invoking the config case command
    result = invoke_cli(
        [
            "config",
            "case",
            "--case-id",
            case_id_tumor_only,
            "--analysis-dir",
            analysis_dir,
            "--fastq-path",
            fastq_dir_tumor_only,
            "--balsamic-cache",
            balsamic_cache,
            "--tumor-sample-name",
            tumor_sample_name,
            "--exome",
        ],
    )
    # THEN config should fail with error message
    assert result.exit_code == 2
    assert "If --exome is provided, --panel-bed must also be provided." in result.output


def test_config_tga_with_exome(
    invoke_cli,
    tumor_sample_name: str,
    analysis_dir: str,
    balsamic_cache: str,
    fastq_dir_tumor_only: str,
    case_id_tumor_only: str,
    panel_bed_file: str,
):
    """Test balsamic config case with GENS arguments for TGA."""

    # GIVEN CLI arguments including optional GENS input-files

    # WHEN invoking the config case command
    result = invoke_cli(
        [
            "config",
            "case",
            "--case-id",
            case_id_tumor_only,
            "--analysis-dir",
            analysis_dir,
            "--fastq-path",
            fastq_dir_tumor_only,
            "--balsamic-cache",
            balsamic_cache,
            "--tumor-sample-name",
            tumor_sample_name,
            "-p",
            panel_bed_file,
            "--exome",
        ],
    )
    # THEN a config should be created and exist
    assert result.exit_code == 0
    assert Path(
        analysis_dir, case_id_tumor_only, f"{case_id_tumor_only}.{FileType.JSON}"
    ).exists()
