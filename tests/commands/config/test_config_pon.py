import json
import graphviz
import logging
from unittest import mock
from pathlib import Path
from BALSAMIC.constants.analysis import PONWorkflow


def test_cnvkit_pon_config(
    invoke_cli,
    analysis_dir: str,
    balsamic_cache: str,
    panel_bed_file: str,
    fastq_dir_pon: str,
    case_id_pon: str,
    sentieon_license: str,
    sentieon_install_dir: str,
):
    """Test balsamic PON config case command for CNVkit."""

    # GIVEN a case ID, fastq files, and an analysis dir

    # WHEN creating a case config
    result = invoke_cli(
        [
            "config",
            "pon",
            "--case-id",
            case_id_pon,
            "--analysis-dir",
            analysis_dir,
            "--fastq-path",
            fastq_dir_pon,
            "-p",
            panel_bed_file,
            "--version",
            "v5",
            "--balsamic-cache",
            balsamic_cache,
            "--pon-workflow",
            PONWorkflow.CNVKIT,
            "--sentieon-install-dir",
            sentieon_install_dir,
            "--sentieon-license",
            sentieon_license,
        ]
    )

    # THEN a config should be created and exist
    assert result.exit_code == 0
    assert Path(analysis_dir, case_id_pon, case_id_pon + "_PON.json").exists()


def test_gens_pon_config(
    invoke_cli,
    analysis_dir: str,
    balsamic_cache: str,
    fastq_dir_gens_pon: str,
    case_id_gens_pon: str,
    gens_hg19_interval_list: str,
    sentieon_license: str,
    sentieon_install_dir: str,
):
    """Test balsamic PON config case command for GENS."""

    # GIVEN a case ID, fastq files, and an analysis dir

    # WHEN creating a config for GENS pon creation workflow
    result = invoke_cli(
        [
            "config",
            "pon",
            "--case-id",
            case_id_gens_pon,
            "--analysis-dir",
            analysis_dir,
            "--fastq-path",
            fastq_dir_gens_pon,
            "--version",
            "v5",
            "--balsamic-cache",
            balsamic_cache,
            "--pon-workflow",
            PONWorkflow.GENS_MALE,
            "--genome-interval",
            gens_hg19_interval_list,
            "--sentieon-install-dir",
            sentieon_install_dir,
            "--sentieon-license",
            sentieon_license,
        ]
    )

    # THEN a config should be created and exist
    assert result.exit_code == 0
    assert Path(analysis_dir, case_id_gens_pon, case_id_gens_pon + "_PON.json").exists()


def test_gens_pon_config(
    invoke_cli,
    analysis_dir: str,
    balsamic_cache: str,
    fastq_dir_gens_pon: str,
    case_id_gens_pon: str,
    sentieon_license: str,
    sentieon_install_dir: str,
):
    """Test detection of missing genome_interval file which is optional but required for GENS."""

    # GIVEN a case ID, fastq files, and an analysis dir

    # WHEN creating a config for GENS pon creation workflow without required genome_interval file
    result = invoke_cli(
        [
            "config",
            "pon",
            "--case-id",
            case_id_gens_pon,
            "--analysis-dir",
            analysis_dir,
            "--fastq-path",
            fastq_dir_gens_pon,
            "--version",
            "v5",
            "--balsamic-cache",
            balsamic_cache,
            "--pon-workflow",
            PONWorkflow.GENS_MALE,
            "--sentieon-install-dir",
            sentieon_install_dir,
            "--sentieon-license",
            sentieon_license,
        ]
    )

    # THEN command should exit with error
    assert (
        "Argument: genome_interval is required for GENS PON creation." in result.output
    )
    assert result.exit_code == 2


def test_cnvkit_pon_config_failed(
    invoke_cli,
    tmp_path: str,
    balsamic_cache: str,
    panel_bed_file: str,
    sentieon_license: str,
    sentieon_install_dir: str,
):
    """Test detection of missing option for a PON config without required arguments.."""
    # GIVEN a case ID, fastq files, and an analysis dir
    test_analysis_dir = tmp_path / "test_analysis_dir"
    test_analysis_dir.mkdir()
    case_id = "sample_pon"

    # WHEN creating config for cnvkit pon creation workflow
    result = invoke_cli(
        [
            "config",
            "pon",
            "--case-id",
            case_id,
            "-p",
            panel_bed_file,
            "--analysis-dir",
            test_analysis_dir,
            "--balsamic-cache",
            balsamic_cache,
            "--sentieon-install-dir",
            sentieon_install_dir,
            "--sentieon-license",
            sentieon_license,
        ]
    )

    # THEN a config should not be created and exit
    assert "Error: Missing option" in result.output
    assert result.exit_code == 2


def test_cnvkit_pon_config_missing_panel(
    invoke_cli,
    tmp_path: str,
    balsamic_cache: str,
    fastq_dir_pon: str,
    sentieon_license: str,
    sentieon_install_dir: str,
):
    """Test detection of missing panel which is optional but required for CNVkit."""

    # GIVEN a case ID, fastq files, and an analysis dir
    test_analysis_dir = tmp_path / "test_analysis_dir"
    test_analysis_dir.mkdir()
    case_id = "sample_pon"

    # WHEN creating config for cnvkit pon creation workflow
    result = invoke_cli(
        [
            "config",
            "pon",
            "--case-id",
            case_id,
            "--analysis-dir",
            test_analysis_dir,
            "--balsamic-cache",
            balsamic_cache,
            "--fastq-path",
            fastq_dir_pon,
            "--pon-workflow",
            PONWorkflow.CNVKIT,
            "--version",
            "v5",
            "--sentieon-install-dir",
            sentieon_install_dir,
            "--sentieon-license",
            sentieon_license,
        ]
    )

    # THEN a config should not be created and exit
    assert "Argument: panel_bed is required for CNVkit PON creation." in result.output
    assert result.exit_code == 2


def test_dag_graph_success_cnvkit_pon(cnvkit_pon_creation_config: str):
    """Test DAG graph building success."""
    # WHEN creating config using standard CLI input and setting Sentieon env vars

    # THEN DAG graph should be created successfully
    assert Path(json.load(open(cnvkit_pon_creation_config))["analysis"]["dag"]).exists()


def test_dag_graph_success_gens_pon(gens_pon_creation_config: str):
    """Test DAG graph building success."""
    # WHEN creating config using standard CLI input and setting Sentieon env vars

    # THEN DAG graph should be created successfully
    assert Path(json.load(open(gens_pon_creation_config))["analysis"]["dag"]).exists()


def test_cnvkit_pon_config_graph_failed(
    invoke_cli,
    analysis_dir: str,
    balsamic_cache: str,
    panel_bed_file: str,
    sentieon_license: str,
    sentieon_install_dir: str,
):
    """Test DAG graph building failure."""

    # GIVEN an analysis config
    case_id = "sample_pon"
    fastq_dir: Path = Path(analysis_dir, case_id, "fastq")
    fastq_dir.mkdir(parents=True, exist_ok=True)

    with mock.patch.object(graphviz, "Source") as mocked:
        mocked.return_value = None
        pon_result = invoke_cli(
            [
                "config",
                "pon",
                "--case-id",
                case_id,
                "--analysis-dir",
                analysis_dir,
                "--fastq-path",
                fastq_dir,
                "-p",
                panel_bed_file,
                "--version",
                "v5",
                "--balsamic-cache",
                balsamic_cache,
                "--pon-workflow",
                PONWorkflow.CNVKIT,
                "--sentieon-install-dir",
                sentieon_install_dir,
                "--sentieon-license",
                sentieon_license,
            ]
        )

    # THEN the graph should not have been built
    assert pon_result.exit_code == 1
