import re
import subprocess
import json
from unittest import mock
from pathlib import Path
from unittest.mock import patch, MagicMock


def test_run_analysis_dragen(invoke_cli, tumor_only_wgs_config):
    # GIVEN a Mock subprocess.run to simulate successful sbatch submission
    mock_result = MagicMock()
    mock_result.stdout = "Submitted batch job 12345"
    mock_result.returncode = 0

    # GIVEN a tumor-only config file
    # WHEN running analysis

    with patch(
        "BALSAMIC.models.sbatchsubmitter.subprocess.run", return_value=mock_result
    ) as mock_run:
        result = invoke_cli(
            ["run", "analysis", "-s", tumor_only_wgs_config, "--dragen"]
        )
        mock_run.assert_called_once()
        assert "sbatch" in mock_run.call_args[0][0]
    assert result.exit_code == 0


def test_run_analysis_tumor_normal_dry_run(invoke_cli, tumor_normal_config):
    # GIVEN a Mock subprocess.run to simulate successful sbatch submission
    mock_result = MagicMock()
    mock_result.stdout = "Submitted batch job 12345"
    mock_result.returncode = 0

    # GIVEN a tumor-normal config file
    # WHEN running analysis

    with patch(
        "BALSAMIC.models.sbatchsubmitter.subprocess.run", return_value=mock_result
    ) as mock_run:
        result = invoke_cli(["run", "analysis", "-s", tumor_normal_config])
        mock_run.assert_called_once()
        assert "sbatch" in mock_run.call_args[0][0]
    # THEN it should run without any error
    assert result.exit_code == 0


def test_run_analysis_tumor_normal_run_interactively(invoke_cli, tumor_normal_config):
    # GIVEN a tumor-normal config file
    # WHEN running analysis
    result = invoke_cli(
        ["run", "analysis", "-s", tumor_normal_config, "--run-interactively"]
    )

    # THEN it should run without error
    assert result.exit_code == 0
    # THEN the start interactive string should be printed
    assert "Starting balsamic workflow interactively" in result.output


def test_run_analysis_tumor_only_dry_run(invoke_cli, tumor_only_config):
    # GIVEN a Mock subprocess.run to simulate successful sbatch submission
    mock_result = MagicMock()
    mock_result.stdout = "Submitted batch job 12345"
    mock_result.returncode = 0

    # GIVEN a tumor-only config file
    # WHEN running analysis

    with patch(
        "BALSAMIC.models.sbatchsubmitter.subprocess.run", return_value=mock_result
    ) as mock_run:
        result = invoke_cli(["run", "analysis", "-s", tumor_only_config])
        mock_run.assert_called_once()
        assert "sbatch" in mock_run.call_args[0][0]
    # THEN it should run without any error
    assert result.exit_code == 0


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
                "--account",
                "development",
            ]
        )
    # THEN it should create a log_dir
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
