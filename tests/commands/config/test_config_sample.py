import json

from unittest import mock

import graphviz

from pathlib import Path

from tests.conftest import MOCKED_OS_ENVIRON


def test_tumor_normal_config(
    invoke_cli,
    case_id_tumor_normal: str,
    tumor_sample_name: str,
    normal_sample_name: str,
    analysis_dir: str,
    fastq_dir_tumor_normal: str,
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
                fastq_dir_tumor_normal,
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
        analysis_dir, case_id_tumor_normal, case_id_tumor_normal + ".json"
    ).exists()


def test_tumor_normal_extrafastq_config(
    invoke_cli,
    case_id_tumor_normal: str,
    tumor_sample_name: str,
    normal_sample_name: str,
    analysis_dir: str,
    fastq_dir_tumor_normal_extrafile: str,
    balsamic_cache: str,
    panel_bed_file: str,
    sentieon_license: str,
    sentieon_install_dir: str,
    caplog,
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
                fastq_dir_tumor_normal_extrafile,
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
    assert (
        "Fastq files found in fastq directory not assigned to any sample" in caplog.text
    )
    assert result.exit_code == 1


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
    assert Path(analysis_dir, case_id_tumor_only, case_id_tumor_only + ".json").exists()


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
    case_id_tumor_only_umi: str,
    tumor_sample_name: str,
    analysis_dir: str,
    fastq_dir_tumor_only_umi: str,
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
            case_id_tumor_only_umi,
            "--analysis-workflow",
            "balsamic-umi",
            "--analysis-dir",
            analysis_dir,
            "--fastq-path",
            fastq_dir_tumor_only_umi,
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
    pon_cnn: str,
    fastq_dir_tumor_only_pon_cnn: str,
    case_id_tumor_only_pon_cnn: str,
):
    """Test balsamic config case with a PON reference."""

    # GIVEN CLI arguments including optional pon reference ".cnn" file

    # WHEN invoking the config case command
    result = invoke_cli(
        [
            "config",
            "case",
            "--case-id",
            case_id_tumor_only_pon_cnn,
            "--analysis-dir",
            analysis_dir,
            "--fastq-path",
            fastq_dir_tumor_only_pon_cnn,
            "-p",
            panel_bed_file,
            "--pon-cnn",
            pon_cnn,
            "--balsamic-cache",
            balsamic_cache,
            "--tumor-sample-name",
            tumor_sample_name,
        ],
    )

    # THEN program exits and checks for filepath
    assert result.exit_code == 0
    assert Path(pon_cnn).exists()


def test_dag_graph_success_tumor_only(tumor_only_config: str):
    """Test DAG graph building success."""
    # WHEN creating config using standard CLI input and setting Sentieon env vars

    # THEN DAG graph should be created successfully
    assert Path(json.load(open(tumor_only_config))["analysis"]["dag"]).exists()


def test_dag_graph_success_tumor_normal(tumor_normal_config: str):
    """Test DAG graph building success."""
    # WHEN creating config using standard CLI input and setting Sentieon env vars

    # THEN DAG graph should be created successfully
    assert Path(json.load(open(tumor_normal_config))["analysis"]["dag"]).exists()


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
