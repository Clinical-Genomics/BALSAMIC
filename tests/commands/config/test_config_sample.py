import os
import json

from unittest import mock

from pathlib import Path
from click.testing import CliRunner
from BALSAMIC.commands.base import cli


def test_dag_graph_success(tumor_normal_wgs_config, tumor_only_config,
                           tumor_normal_config, tumor_only_wgs_config,
                           tumor_only_umi_config):
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


def test_tumor_only_config_bad_filename(
        tmp_path_factory,
        analysis_dir,
        panel_bed_file,
        balsamic_cache
):

    # GIVEN existing fastq file with wrong naming convention
    faulty_fastq_dir = tmp_path_factory.mktemp("error_fastq")
    Path(faulty_fastq_dir / "error.fastq.gz").touch()

    case_id = "faulty_tumor"
    tumor = Path(faulty_fastq_dir / "error.fastq.gz").as_posix()

    # Invoke CLI command using file as argument
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "config",
            "case",
            "-t",
            tumor,
            "-p",
            panel_bed_file,
            "--case-id",
            case_id,
            "--analysis-dir",
            analysis_dir,
"--balsamic-cache",
balsamic_cache,
        ],
    )

    # THEN run should abort
    assert result.exit_code == 1


def test_run_without_permissions(
        no_write_perm_path,
        sample_fastq,
        panel_bed_file,
        balsamic_cache,
):
    # GIVEN CLI arguments including an analysis_dir without write permissions
    case_id = "sample_tumor_only"
    tumor = sample_fastq["tumor"]

    runner = CliRunner()
    result = runner.invoke(
        cli,
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
        sample_fastq, analysis_dir, balsamic_cache,
        panel_bed_file):

    # GIVEN CLI arguments including a background variant file
    case_id = "sample_umi_tumor_only"
    tumor = sample_fastq["tumor"]
    background_file = "tests/test_data/references/panel/background_variants.txt"
    background_variant_file = background_file

    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "config", "case", "-p", panel_bed_file, "-t", tumor, "--case-id",
            case_id, "--analysis-dir", analysis_dir, 
"--balsamic-cache",
balsamic_cache,
            "--background-variants", background_variant_file
        ],
    )
    # THEN program exits and checks for filepath
    assert result.exit_code == 0
    assert Path(background_variant_file).exists()
