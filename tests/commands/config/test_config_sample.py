import os
import json

from unittest import mock

import graphviz

from pathlib import Path


def test_tumor_normal_config(invoke_cli, sample_fastq, tmp_path,
                        balsamic_cache, panel_bed_file,
                        sentieon_license, sentieon_install_dir):
    # GIVEN a case ID, fastq files, and an analysis dir
    test_analysis_dir = tmp_path / "test_analysis_dir"
    test_analysis_dir.mkdir()
    case_id = "sample_tumor_normal"
    tumor = sample_fastq["tumor"]
    normal = sample_fastq["normal"]

    # WHEN creating a case analysis
    with mock.patch.dict(
            "os.environ", {
                'SENTIEON_LICENSE': sentieon_license,
                'SENTIEON_INSTALL_DIR': sentieon_install_dir
            }):
        result = invoke_cli(
            [
                "config",
                "case",
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
    assert Path(test_analysis_dir, case_id, case_id + ".json").exists()


def test_tumor_only_config(invoke_cli, sample_fastq, tmp_path,
                        balsamic_cache, panel_bed_file,
                        sentieon_license, sentieon_install_dir):
    # GIVEN a case ID, fastq files, and an analysis dir
    test_analysis_dir = tmp_path / "test_analysis_dir"
    test_analysis_dir.mkdir()
    case_id = "sample_tumor_only"
    tumor = sample_fastq["tumor"]

    # WHEN creating a case analysis
    with mock.patch.dict(
            "os.environ", {
                'SENTIEON_LICENSE': sentieon_license,
                'SENTIEON_INSTALL_DIR': sentieon_install_dir
            }):
        result = invoke_cli(
            [
                "config",
                "case",
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
    assert Path(test_analysis_dir, case_id, case_id + ".json").exists()

def test_dag_graph_success(tumor_normal_wgs_config, tumor_only_config,
                           tumor_normal_config, tumor_only_wgs_config,
                           tumor_only_umi_config, pon_config):
    # WHEN creating config using standard CLI input and setting Sentieon env vars
    # THEN DAG graph should be created successfully
    assert Path(json.load(
        open(tumor_normal_config))["analysis"]["dag"]).exists()
    assert Path(json.load(open(tumor_only_config))["analysis"]["dag"]).exists()
    assert Path(json.load(
        open(tumor_only_wgs_config))["analysis"]["dag"]).exists()
    assert Path(json.load(
        open(tumor_normal_wgs_config))["analysis"]["dag"]).exists()
    assert Path(json.load(
        open(tumor_only_umi_config))["analysis"]["dag"]).exists()
    assert Path(json.load(
        open(pon_config))["analysis"]["dag"]).exists()


def test_config_bad_filename(
        invoke_cli,
        tmp_path_factory,
        analysis_dir,
        panel_bed_file,
        balsamic_cache,
):
    # GIVEN existing fastq file with wrong naming convention
    faulty_fastq_dir = tmp_path_factory.mktemp("error_fastq")
    Path(faulty_fastq_dir / "error.fastq.gz").touch()

    case_id1 = "faulty_tumor"
    case_id2 = "faulty_pon_normal"
    tumor = Path(faulty_fastq_dir / "error.fastq.gz").as_posix()
    normal = Path(faulty_fastq_dir / "error.fastq.gz").as_posix()


    # Invoke CLI command using file as argument
    case_result = invoke_cli(
        [
            "config",
            "case",
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
    assert case_result.exit_code == 1
    assert pon_result.exit_code == 1


def test_run_without_permissions(
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
            "case",
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


def test_tumor_only_umi_config_background_file(
        invoke_cli,
        sample_fastq, analysis_dir, balsamic_cache, 
        panel_bed_file):

    # GIVEN CLI arguments including a background variant file
    case_id = "sample_umi_tumor_only"
    tumor = sample_fastq["tumor"]
    background_file = "tests/test_data/references/panel/background_variants.txt"
    background_variant_file = background_file

    result = invoke_cli(
        [
            "config", "case", "-p", panel_bed_file, "-t", tumor, "--case-id",
            case_id, "--analysis-dir", analysis_dir, 
            "--background-variants", background_variant_file, "--balsamic-cache",
            balsamic_cache,
        ],
    )
    # THEN program exits and checks for filepath
    assert result.exit_code == 0
    assert Path(background_variant_file).exists()


def test_config_graph_failed(invoke_cli, sample_fastq,
                                  analysis_dir, balsamic_cache,
                                  panel_bed_file):
    # GIVEN an analysis config 
    case_id = "sample_tumor_only"
    pon_case_id = "sample_pon"
    tumor = sample_fastq["tumor"]
    normal = sample_fastq["normal"]
    normal1 = sample_fastq["normal1"]
    normal2 = sample_fastq["normal2"]

    with mock.patch.object(graphviz, 'Source') as mocked:
        mocked.return_value = None
        case_result = invoke_cli(
            [
                "config",
                "case",
                "-p",
                panel_bed_file,
                "-t",
                tumor,
                "--case-id",
                case_id,
                "--analysis-dir",
                analysis_dir,
                "--balsamic-cache",
                balsamic_cache
            ],
        )
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

    assert case_result.exit_code == 1
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