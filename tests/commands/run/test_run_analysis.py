import re
import subprocess
import json
from unittest import mock
from pathlib import Path
from pydantic import ValidationError



def test_run_analysis_dragen(invoke_cli, tumor_only_wgs_config):
    # GIVEN a WGS config file
    # WHEN running analysis
    result = invoke_cli(["run", "analysis", "-s", tumor_only_wgs_config, "--dragen"])

    # THEN it should run without any error
    assert result.exit_code == 0


def test_run_analysis_disable_variant_caller(invoke_cli, tumor_only_config):
    # GIVEN a tumor-only config file and variant caller to disable
    disabled_varcaller = "mutect"

    # WHEN running analysis
    result = invoke_cli(
        [
            "run",
            "analysis",
            "-s",
            tumor_only_config,
            "--disable-variant-caller",
            disabled_varcaller,
        ]
    )

    # THEN it should run without any error
    assert result.exit_code == 0
    assert disabled_varcaller not in result.output


def test_run_analysis_tumor_normal_dry_run(invoke_cli, tumor_normal_config):
    # GIVEN a tumor-normal config file
    # WHEN running analysis
    result = invoke_cli(["run", "analysis", "-s", tumor_normal_config])

    # THEN it should run without any error
    assert result.exit_code == 0


def test_run_analysis_tumor_only_dry_run(
    invoke_cli, tumor_only_config, tumor_normal_config
):
    # GIVEN a tumor-only config file
    # WHEN running analysis
    result = invoke_cli(
        ["run", "analysis", "-s", tumor_only_config, "--mail-type", "FAIL"]
    )

    # THEN it should run without any error
    assert result.exit_code == 0

def test_validate_fastq_input_illegal_pairnames(invoke_cli, sample_config_illegal_fastq_pairnames):
    """Tests ability of balsamic model to detect non-existent fastq files assigned to sample config."""

    # GIVEN a tumor normal config where fastq_dir is empty

    # WHEN instantiating the balsamic model
    # THEN the following error should be found

    result = invoke_cli(
        ["run", "analysis", "-s", sample_config_illegal_fastq_pairnames]
    )
    assert result.exit_code == 1
    exception = result.exception
    assert isinstance(exception, ValidationError)
    error_message = str(exception)
    assert (
            "Fastq pair does not have names of equal length"
            in error_message
    )

def test_validate_fastq_input_missingfiles(invoke_cli, sample_config_new_data):
    """Tests ability of balsamic model to detect non-existent fastq files assigned to sample config."""

    # GIVEN a tumor normal config where fastq_dir is empty

    # WHEN instantiating the balsamic model
    # THEN the following error should be found

    result = invoke_cli(
        ["run", "analysis", "-s", sample_config_new_data]
    )
    assert result.exit_code == 1
    exception = result.exception
    assert isinstance(exception, ValidationError)
    error_message = str(exception)
    assert (
            "Fastqs in fastq-dir not assigned to sample config"
            in error_message
    )


def test_run_analysis_click_abort(invoke_cli, tumor_only_config, tumor_normal_config):
    # GIVEN a tumor-only config file
    # WHEN running analysis
    result = invoke_cli(["run", "analysis", "-r", "-s", tumor_only_config])

    # THEN it should abort with error
    assert result.exit_code == 1


def test_run_analysis_create_dir(invoke_cli, tumor_only_config):
    # GIVEN a tumor-only config file
    # WHEN running analysis

    with open(tumor_only_config) as fh:
        tumor_config = json.load(fh)
    log_dir = tumor_config["analysis"]["log"]
    Path(log_dir).mkdir(exist_ok=True)
    Path(log_dir, "logfile.log").touch(exist_ok=True)

    with mock.patch.object(subprocess, "run") as mocked:
        mocked.return_value.stdout = 1
        invoke_cli(
            [
                "run",
                "analysis",
                "-s",
                tumor_only_config,
                "-r",
                "--benchmark",
                "--account",
                "development",
            ]
        )
        # THEN it should abort with error
        assert Path(re.sub("/$", ".1/", log_dir)).exists()


def test_run_analysis_ponpath(invoke_cli, tumor_only_pon_config):
    # GIVEN a tumor-only with pon file in the config file
    # WHEN running analysis

    with open(tumor_only_pon_config) as fh:
        sample_config = json.load(fh)

    bind_path = ["/path_to_dummy/ash/"]
    pon_fl = sample_config["panel"].get("pon_cnn")
    pon_path = Path(pon_fl).resolve()

    if "pon_cnn" in sample_config["panel"]:
        bind_path.append(str(pon_path))

    # THEN it checks for existence of paths
    assert pon_path.exists()
    assert str(pon_path) in bind_path
