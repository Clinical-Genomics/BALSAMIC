import os
import json
import graphviz
import logging
from unittest import mock
from pathlib import Path

from BALSAMIC.utils.cli import create_pon_fastq_symlink

qc_json = '_QC.json'

def test_qc_normal_config(
    invoke_cli,
    sample_fastq,
    tmp_path,
    balsamic_cache,
    panel_bed_file,
):
    # GIVEN a case ID, fastq files, and an analysis dir
    test_analysis_dir = tmp_path / "test_analysis_dir"
    test_analysis_dir.mkdir()
    case_id = "sample_tumor_normal"
    tumor = sample_fastq["tumor"]
    normal = sample_fastq["normal"]

    # WHEN creating a case analysis
    with mock.patch.dict(
        "os.environ",
    ):
        result = invoke_cli(
            [
                "config",
                "qc",
                "-p",
                panel_bed_file,
                "-t",
                tumor,
                "-n",
                normal,
                "--case-id",
                case_id,
                "--analysis-dir",
                test_analysis_dir,
                "--balsamic-cache",
                balsamic_cache,
                "--tumor-sample-name",
                "ACC1",
                "--normal-sample-name",
                "ACC2",
            ],
        )

    # THEN a config should be created and exist
    assert result.exit_code == 0
    assert Path(test_analysis_dir, case_id, case_id + qc_json).exists()
    # load json file and check if dag exists
    qc_config = json.load(open(Path(test_analysis_dir, case_id, case_id + qc_json)))
    # assert if config json dag file is created
    assert Path(qc_config["analysis"]["dag"]).exists()


def test_qc_tumor_only_config(
    invoke_cli,
    sample_fastq,
    tmp_path,
    balsamic_cache,
    panel_bed_file,
    sentieon_license,
    sentieon_install_dir,
):
    # GIVEN a case ID, fastq files, and an analysis dir
    test_analysis_dir = tmp_path / "test_analysis_dir"
    test_analysis_dir.mkdir()
    case_id = "sample_tumor_only"
    tumor = sample_fastq["tumor"]

    # WHEN creating a case analysis
    with mock.patch.dict(
        "os.environ",
        {
            "SENTIEON_LICENSE": sentieon_license,
            "SENTIEON_INSTALL_DIR": sentieon_install_dir,
        },
    ):
        result = invoke_cli(
            [
                "config",
                "qc",
                "-p",
                panel_bed_file,
                "-t",
                tumor,
                "--case-id",
                case_id,
                "--analysis-dir",
                test_analysis_dir,
                "--balsamic-cache",
                balsamic_cache,
                "--tumor-sample-name",
                "ACC1",
            ],
        )

    # THEN a config should be created and exist
    assert result.exit_code == 0
    assert Path(test_analysis_dir, case_id, case_id + qc_json).exists()
    # load json file and check if dag exists
    qc_config = json.load(open(Path(test_analysis_dir, case_id, case_id + qc_json)))
    # assert if config json dag file is created
    assert Path(qc_config["analysis"]["dag"]).exists()


def test_qc_config_bad_filename(
    invoke_cli,
    tmp_path_factory,
    analysis_dir,
    panel_bed_file,
    balsamic_cache,
):
    # GIVEN existing fastq file with wrong naming convention
    faulty_fastq_dir = tmp_path_factory.mktemp("error_fastq")
    fastq_file_name_tumor = "tumor_error.fastq.gz"
    Path(faulty_fastq_dir / fastq_file_name_tumor).touch()

    case_id1 = "faulty_tumor"
    tumor = Path(faulty_fastq_dir / fastq_file_name_tumor).as_posix()

    # Invoke CLI command using file as argument
    case_result = invoke_cli(
        [
            "config",
            "qc",
            "-t",
            tumor,
            "-p",
            panel_bed_file,
            "--case-id",
            case_id1,
            "--analysis-dir",
            analysis_dir,
            "--balsamic-cache",
            balsamic_cache,
        ],
    )

    # THEN run should abort
    assert case_result.exit_code == 1


def test_qc_run_without_permissions(
    invoke_cli,
    no_write_perm_path,
    sample_fastq,
    panel_bed_file,
    balsamic_cache,
):
    # GIVEN CLI arguments including an analysis_dir without write permissions
    case_id = "sample_tumor_only"
    tumor = sample_fastq["tumor"]

    result = invoke_cli(
        [
            "config",
            "qc",
            "-p",
            panel_bed_file,
            "-t",
            tumor,
            "--case-id",
            case_id,
            "--analysis-dir",
            no_write_perm_path,
            "--balsamic-cache",
            balsamic_cache,
        ],
    )
    # THEN program exits before completion
    assert result.exit_code == 1


def test_qc_config_failed(invoke_cli, tmp_path, balsamic_cache, panel_bed_file):
    # GIVEN a case ID, fastq files, and an analysis dir
    test_analysis_dir = tmp_path / "test_analysis_dir"
    test_analysis_dir.mkdir()
    case_id = "sample_pon"

    # WHEN creating a case analysis
    result = invoke_cli(
        [
            "config",
            "qc",
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


def test_config_qc_graph_failed(
    invoke_cli, sample_fastq, analysis_dir, balsamic_cache, panel_bed_file
):
    # GIVEN an analysis config
    case_id = "sample_tumor_only"
    tumor = sample_fastq["tumor"]

    with mock.patch.object(graphviz, "Source") as mocked:
        mocked.return_value = None
        case_result = invoke_cli(
            [
                "config",
                "qc",
                "-p",
                panel_bed_file,
                "-t",
                tumor,
                "--case-id",
                case_id,
                "--analysis-dir",
                analysis_dir,
                "--balsamic-cache",
                balsamic_cache,
            ],
        )

    assert case_result.exit_code == 1
