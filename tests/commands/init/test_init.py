"""Test Balsamic init command."""
from functools import partial
from pathlib import Path
from unittest.mock import patch, MagicMock
from click.testing import Result
from graphviz import Source

from BALSAMIC import __version__ as balsamic_version
from BALSAMIC.constants.analysis import RunMode
from BALSAMIC.constants.cache import GenomeVersion
from BALSAMIC.constants.cluster import ClusterAccount
from BALSAMIC.constants.constants import EXIT_SUCCESS, EXIT_FAIL
from BALSAMIC.utils.exc import BalsamicError
from BALSAMIC.models.sbatchsubmitter import SbatchSubmitter


def test_init_hg(
    invoke_cli: partial,
    tmp_path: Path,
    cosmic_key: str,
    config_json: str,
    reference_graph: str,
):
    """Test Balsamic init command."""

    # GIVEN a temporary output directory and a COSMIC key

    # WHEN invoking the init command
    result: Result = invoke_cli(
        [
            "init",
            "--out-dir",
            tmp_path.as_posix(),
            "--genome-version",
            GenomeVersion.HG19,
            "--cosmic-key",
            cosmic_key,
            "--run-interactively",
        ]
    )

    # THEN the human reference generation workflow should have successfully started
    assert Path(tmp_path, balsamic_version, GenomeVersion.HG19, config_json).exists()
    assert Path(
        tmp_path, balsamic_version, GenomeVersion.HG19, reference_graph
    ).exists()
    assert result.exit_code == EXIT_SUCCESS


def test_init_hg_no_cosmic_key(invoke_cli: partial, tmp_path: Path, cosmic_key: str):
    """Test Balsamic init command when a COSMIC key is not provided."""

    # GIVEN a temporary output directory and a COSMIC key

    # WHEN invoking the init command
    result: Result = invoke_cli(
        [
            "init",
            "--out-dir",
            tmp_path.as_posix(),
            "--genome-version",
            GenomeVersion.HG19,
        ]
    )

    # THEN an exception should have been raised
    assert (
        f"No COSMIC authentication key specified. It is required when using {GenomeVersion.HG19} reference"
        in result.output
    )
    assert result.exit_code == EXIT_FAIL


def test_init_hg_run_analysis_no_account(
    invoke_cli: partial, tmp_path: Path, cosmic_key: str
):
    """Test Balsamic init command when actually running the analysis without specifying a cluster account."""

    # GIVEN a temporary output directory and a COSMIC key

    # WHEN invoking the init command
    result: Result = invoke_cli(
        [
            "init",
            "--out-dir",
            tmp_path.as_posix(),
            "--genome-version",
            GenomeVersion.HG19,
            "--cosmic-key",
            cosmic_key,
            "--run-mode",
            RunMode.CLUSTER,
            "--run-analysis",
        ]
    )

    # THEN an exception should have been raised
    assert "A cluster account is required for cluster run mode" in result.output
    assert result.exit_code == EXIT_FAIL


def test_init_hg_submit_succeeds(
    invoke_cli: partial,
    tmp_path: Path,
    cosmic_key: str,
):
    # Mock sbatch returning a job id
    mock_sbatch = MagicMock()
    mock_sbatch.stdout = "Submitted batch job 12345"
    mock_sbatch.returncode = 0

    with (
        # IMPORTANT: skip the pre-step graph generation that shells out
        patch("BALSAMIC.commands.init.base.generate_graph") as mock_gen_graph,
        # Patch the *submitter* subprocess.run (this is where 'sbatch' is called)
        patch(
            "BALSAMIC.models.sbatchsubmitter.subprocess.run", return_value=mock_sbatch
        ) as mock_submit_run,
        # Optional: if you want to assert write_job_id_yaml call
        patch(
            "BALSAMIC.commands.init.base.SbatchSubmitter.write_job_id_yaml"
        ) as mock_write_yaml,
    ):
        result: Result = invoke_cli(
            [
                "init",
                "--out-dir",
                tmp_path.as_posix(),
                "--genome-version",
                GenomeVersion.HG19,
                "--cosmic-key",
                cosmic_key,
                # Do NOT pass any flag that turns on interactive mode
            ]
        )

    # `generate_graph` was called (but we stubbed it out)
    mock_gen_graph.assert_called_once()

    # The submitter’s subprocess.run should have been called once with ["sbatch", ...]
    mock_submit_run.assert_called_once()
    cmd = mock_submit_run.call_args[0][0]
    assert isinstance(cmd, list) and cmd and cmd[0] == "sbatch"

    # When job id is present, we expect to write it
    mock_write_yaml.assert_called_once_with("12345")

    assert result.exit_code == EXIT_SUCCESS


def test_init_hg_submit_no_jobid_logs_warning(
    invoke_cli: partial,
    tmp_path: Path,
    cosmic_key: str,
):
    # Simulate sbatch output that does NOT contain a job id
    mock_sbatch = MagicMock()
    mock_sbatch.stdout = "weird output without job id"
    mock_sbatch.returncode = 0

    with (
        patch("BALSAMIC.commands.init.base.generate_graph"),
        patch(
            "BALSAMIC.models.sbatchsubmitter.subprocess.run", return_value=mock_sbatch
        ),
        patch(
            "BALSAMIC.commands.init.base.SbatchSubmitter.write_job_id_yaml"
        ) as mock_write_yaml,
        patch("BALSAMIC.commands.init.base.LOG") as mock_log,
    ):
        result = invoke_cli(
            [
                "init",
                "--out-dir",
                tmp_path.as_posix(),
                "--genome-version",
                GenomeVersion.HG19,
                "--cosmic-key",
                cosmic_key,
            ]
        )

    mock_write_yaml.assert_not_called()
    mock_log.warning.assert_any_call("Could not retrieve job id from SLURM.")
    assert result.exit_code == EXIT_SUCCESS
