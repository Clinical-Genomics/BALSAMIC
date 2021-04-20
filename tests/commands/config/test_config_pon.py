import os
import json

from unittest import mock

import graphviz

from pathlib import Path


def test_config_pon_bad_filename(
        invoke_cli,
        tmp_path_factory,
        analysis_dir,
        panel_bed_file,
        balsamic_cache,
):
    # GIVEN existing fastq file with wrong naming convention
    faulty_fastq_dir = tmp_path_factory.mktemp("error_fastq")
    fastq_file_name_normal = "normal_error.fastq.gz"
    Path(faulty_fastq_dir / fastq_file_name_normal).touch()
    case_id2 = "faulty_pon_normal"
    normal = Path(faulty_fastq_dir / fastq_file_name_normal).as_posix()

    # Invoke CLI command using file as argument
    pon_result = invoke_cli(
        [
            "config",
            "pon",
            "-n",
            normal,
            "-p",
            panel_bed_file,
            "--case-id",
            case_id2,
            "--analysis-dir",
            analysis_dir,
            "--balsamic-cache",
            balsamic_cache,
        ],
    )

    # THEN run should abort
    assert pon_result.exit_code == 1




def test_config_pon_graph_failed(invoke_cli, sample_fastq,
                                 analysis_dir, balsamic_cache,
                                 panel_bed_file):
    # GIVEN an analysis config
    pon_case_id = "sample_pon"
    normal = sample_fastq["normal"]
    normal1 = sample_fastq["normal1"]
    normal2 = sample_fastq["normal2"]

    with mock.patch.object(graphviz, 'Source') as mocked:
        mocked.return_value = None
        pon_result = invoke_cli(
            [
                "config",
                "pon",
                "-p",
                panel_bed_file,
                "-n",
                normal,
                "-n",
                normal1,
                "-n",
                normal2,
                "--case-id",
                pon_case_id,
                "--analysis-dir",
                analysis_dir,
                "--balsamic-cache",
                balsamic_cache
            ],
        )

    assert pon_result.exit_code == 1


def test_pon_config(invoke_cli, sample_fastq, tmp_path,
                    balsamic_cache, panel_bed_file):
    # GIVEN a case ID, fastq files, and an analysis dir
    test_analysis_dir = tmp_path / "test_analysis_dir"
    test_analysis_dir.mkdir()
    case_id = "sample_pon"
    normal = sample_fastq["normal"]
    normal1 = sample_fastq["normal1"]
    normal2 = sample_fastq["normal2"]

    # WHEN creating a case config
    result = invoke_cli(
        [
            "config",
            "pon",
            "--case-id",
            case_id,
            "-p",
            panel_bed_file,
            "-n",
            normal,
            "-n",
            normal1,
            "-n",
            normal2,
            "--analysis-dir",
            test_analysis_dir,
            "--balsamic-cache",
            balsamic_cache
        ]
    )

    # THEN a config should be created and exist
    assert result.exit_code == 0
    assert Path(test_analysis_dir, case_id, case_id + "_PON.json").exists()


def test_pon_config_failed(invoke_cli, tmp_path,
                           balsamic_cache, panel_bed_file):
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
            balsamic_cache
        ]
    )

    # THEN a config should be created and exist
    assert 'Error: Missing option' in result.output
    assert result.exit_code == 2


def test_pon_dag_graph_success(pon_config):
    # WHEN creating config using standard CLI input
    # THEN DAG graph should be created successfully
    assert Path(json.load(
        open(pon_config))["analysis"]["dag"]).exists()
