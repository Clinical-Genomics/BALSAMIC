import json
import graphviz
import logging
from unittest import mock
from pathlib import Path

from BALSAMIC.utils.cli import create_pon_fastq_symlink


def test_pon_config(
    invoke_cli,
    analysis_dir: str,
    balsamic_cache: str,
    panel_bed_file: str,
    pon_fastq_path: str,
):
    """Test balsamic PON config case command."""

    # GIVEN a case ID, fastq files, and an analysis dir
    case_id = "sample_pon"
    fastq_dir: Path = Path(analysis_dir, case_id, "fastq")
    fastq_dir.mkdir(parents=True, exist_ok=True)

    # WHEN creating a case config
    result = invoke_cli(
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
        ]
    )

    # THEN a config should be created and exist
    assert result.exit_code == 0
    assert Path(analysis_dir, case_id, case_id + "_PON.json").exists()


def test_pon_config_failed(invoke_cli, tmp_path, balsamic_cache, panel_bed_file):
    # GIVEN a case ID, fastq files, and an analysis dir
    test_analysis_dir = tmp_path / "test_analysis_dir"
    test_analysis_dir.mkdir()
    case_id = "sample_pon"

    # WHEN creating a case analysis
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
        ]
    )

    # THEN a config should not be created and exit
    assert "Error: Missing option" in result.output
    assert result.exit_code == 2


def test_create_pon_fastq_symlink(tmp_path_factory, caplog):
    pon_symlink_from = tmp_path_factory.mktemp("pon_symlink_from")
    pon_symlink_to = tmp_path_factory.mktemp("pon_symlink_to")
    files = ["normal1_R_1.fastq.gz", "normal1_R_2.fastq.gz"]
    ponfiles = [Path(pon_symlink_from, x) for x in files]
    for ponfile in ponfiles:
        ponfile.touch()
    with caplog.at_level(logging.INFO):
        create_pon_fastq_symlink(pon_symlink_from, pon_symlink_to)
        # THEN destination should have files
        assert len(list(Path(pon_symlink_to).rglob("*.fastq.gz"))) == 2
        # THEN exception triggers log message containing "skipping"
        assert "skipping" not in caplog.text


def test_pon_config_graph(
    invoke_cli,
    tumor_only_pon_config: str,
):
    """Test DAG graph building success."""

    # WHEN creating config using standard CLI input and setting Sentieon env vars

    # THEN DAG graph should be created successfully
    assert Path(json.load(open(tumor_only_pon_config))["analysis"]["dag"]).exists()


def test_pon_config_graph_failed(
    invoke_cli,
    analysis_dir: str,
    balsamic_cache: str,
    panel_bed_file: str,
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
            ]
        )

    # THEN the graph should not have been built
    assert pon_result.exit_code == 1
