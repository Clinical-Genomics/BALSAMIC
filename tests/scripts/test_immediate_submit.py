"""Test immediate submit script."""
import logging
import subprocess
from pathlib import Path
from typing import Any, Dict
from unittest import mock

from _pytest.logging import LogCaptureFixture
from click.testing import CliRunner, Result
from snakemake import utils

from BALSAMIC.assets.scripts.immediate_submit import immediate_submit
from BALSAMIC.constants.constants import EXIT_SUCCESS


def test_immediate_submit(
    job_id: str,
    job_properties: Dict[str, Any],
    scheduler_data: Dict[str, Any],
    session_tmp_path: Path,
    caplog: LogCaptureFixture,
    cli_runner: CliRunner,
):
    """Test immediate submit script execution."""
    caplog.set_level(logging.INFO)

    # GIVEN some scheduler data and a CLI runner with the expected standard output
    stdout: str = f"Submitted batch job {job_id}"

    # WHEN calling the immediate submit script
    with mock.patch.object(
        utils, "read_job_properties", return_value=job_properties
    ), mock.patch.object(
        subprocess,
        "run",
        return_value=subprocess.CompletedProcess(
            args="", returncode=EXIT_SUCCESS, stdout=stdout, stderr=""
        ),
    ):
        result: Result = cli_runner.invoke(
            immediate_submit,
            [
                scheduler_data["case_id"],
                "dependencies",
                scheduler_data["job_script"],
                "--account",
                "development",
                "--log-dir",
                scheduler_data["log_dir"],
                "--mail-user",
                "balsamic@scilifelab.se",
                "--profile",
                "slurm",
                "--qos",
                "high",
                "--script-dir",
                session_tmp_path.as_posix(),
            ],
        )

    # THEN the command should succeed and the log file created
    assert result.exit_code == 0
    assert f"Submitted job with ID: {job_id}" in caplog.text
    assert Path(
        scheduler_data["log_dir"], f"{scheduler_data['case_id']}.sacct"
    ).is_file()
