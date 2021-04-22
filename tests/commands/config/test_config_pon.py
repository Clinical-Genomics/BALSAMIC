import os
import json

from unittest import mock

import graphviz

from pathlib import Path


def test_pon_config(invoke_cli, sample_fastq, tmp_path, balsamic_cache,
                    panel_bed_file, pon_fastq_path):
    # GIVEN a case ID, fastq files, and an analysis dir
    case_id = "sample_pon"
    test_analysis_dir = tmp_path / "test_analysis_dir"
    test_analysis_dir.mkdir()

    # WHEN creating a case config
    result = invoke_cli([
        "config", "pon", "--case-id", case_id, "-p", panel_bed_file,
        "--analysis-dir", test_analysis_dir, "--fastq-path", pon_fastq_path,
        "--balsamic-cache", balsamic_cache
    ])

    # THEN a config should be created and exist
    assert result.exit_code == 0
    assert Path(test_analysis_dir, case_id, case_id + "_PON.json").exists()


def test_pon_config_failed(invoke_cli, tmp_path, balsamic_cache,
                           panel_bed_file):
    # GIVEN a case ID, fastq files, and an analysis dir
    test_analysis_dir = tmp_path / "test_analysis_dir"
    test_analysis_dir.mkdir()
    case_id = "sample_pon"

    # WHEN creating a case analysis
    result = invoke_cli([
        "config", "pon", "--case-id", case_id, "-p", panel_bed_file,
        "--analysis-dir", test_analysis_dir, "--balsamic-cache", balsamic_cache
    ])

    # THEN a config should be created and exist
    assert 'Error: Missing option' in result.output
    assert result.exit_code == 2


def test_config_pon_graph_failed(invoke_cli, sample_fastq, analysis_dir,
                                 balsamic_cache, pon_fastq_path,
                                 panel_bed_file):
    # GIVEN an analysis config
    pon_case_id = "sample_pon"

    with mock.patch.object(graphviz, 'Source') as mocked:
        mocked.return_value = None
        pon_result = invoke_cli([
            "config", "pon", "-p", panel_bed_file, "--case-id", pon_case_id,
            "--analysis-dir", analysis_dir, "--balsamic-cache", balsamic_cache,
            "--fastq-path", pon_fastq_path
        ], )

    assert pon_result.exit_code == 1
