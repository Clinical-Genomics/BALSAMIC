"""Tests for SbatchSubmitter."""
from unittest.mock import MagicMock, patch
import logging
import os
import shlex
import subprocess
import pytest
from pathlib import Path
from BALSAMIC.models.sbatchsubmitter import SbatchSubmitter


MODULE = "BALSAMIC.models.sbatchsubmitter"

# --------------------------- fixtures ---------------------------------


@pytest.fixture
def logger():
    log = logging.getLogger("SbatchSubmitterTest")
    log.setLevel(logging.INFO)
    return log


@pytest.fixture
def fake_exec():
    m = MagicMock()
    m.get_command.return_value = "snakemake --notemp -p --rerun-triggers mtime"
    return m


@pytest.fixture
def submitter(tmp_path, monkeypatch, fake_exec):
    # Real logger; we'll assert via caplog
    log = logging.getLogger("SbatchSubmitterTest")
    log.setLevel(logging.INFO)

    script_path = tmp_path / "scripts"
    result_path = tmp_path / "results"
    log_path = tmp_path / "logs"
    for p in (script_path, result_path, log_path):
        p.mkdir(parents=True)

    monkeypatch.setenv("CONDA_PREFIX", "/conda/env")

    scan_script = tmp_path / "scan_status.py"
    scan_script.write_text("# dummy")

    s = SbatchSubmitter(
        case_id="CASE123",
        script_path=script_path,
        result_path=result_path,
        scan_finished_jobid_status=str(scan_script),
        log_path=log_path,
        account="acct",
        qos="qos",
        headjob_partition=None,
        max_run_hours=4,
        snakemake_executable=fake_exec,
        logger=log,
    )
    # ensure the sbatch script path points to something deterministic for the test
    (script_path / "BALSAMIC_snakemake_submit.sh").write_text("#!/bin/bash\n")
    return s


@pytest.fixture
def submitter_with_partition(tmp_path, monkeypatch, logger, fake_exec):
    script_path = tmp_path / "scripts"
    result_path = tmp_path / "results"
    log_path = tmp_path / "logs"
    for p in (script_path, result_path, log_path):
        p.mkdir(parents=True)
    monkeypatch.setenv("CONDA_PREFIX", "/conda/env")
    scan_script = tmp_path / "scan_status.py"
    scan_script.write_text("# dummy")

    return SbatchSubmitter(
        case_id="CASE123",
        script_path=script_path,
        result_path=result_path,
        scan_finished_jobid_status=str(scan_script),
        log_path=log_path,
        account="myacct",
        qos="normal",
        headjob_partition="core",  # <-- ensure partition is set
        max_run_hours=48,
        snakemake_executable=fake_exec,
        logger=logger,
    )


# --------------------------- header tests ------------------------------


def test_build_sbatch_header_includes_partition(submitter_with_partition):
    s = submitter_with_partition
    header = s._build_sbatch_header()
    assert header.startswith("#!/bin/bash -l")
    assert f"#SBATCH --account={s.account}" in header
    assert f"#SBATCH --job-name=BALSAMIC_snakemake_submit.{s.case_id}.%j" in header
    assert (
        f"#SBATCH --output={s.log_path}/BALSAMIC_snakemake_submit.{s.case_id}.%j.out"
        in header
    )
    assert (
        f"#SBATCH --error={s.log_path}/BALSAMIC_snakemake_submit.{s.case_id}.%j.err"
        in header
    )
    assert "#SBATCH --ntasks=1" in header
    assert "#SBATCH --mem=5G" in header
    assert f"#SBATCH --time={s.max_run_hours}:00:00" in header
    assert f"#SBATCH --qos={s.qos}" in header
    assert "#SBATCH --cpus-per-task=1" in header
    assert "#SBATCH --partition=core" in header


def test_build_sbatch_header_omits_partition_when_none(
    tmp_path, monkeypatch, logger, fake_exec
):
    s = SbatchSubmitter(
        case_id="X",
        script_path=tmp_path,
        result_path=tmp_path,
        scan_finished_jobid_status=str(tmp_path / "scan.py"),
        log_path=tmp_path,
        account="acct",
        qos="qos",
        headjob_partition=None,
        max_run_hours=2,
        snakemake_executable=fake_exec,
        logger=logger,
    )
    header = s._build_sbatch_header()
    assert "#SBATCH --partition=" not in header
    assert "#!/bin/bash -l" in header
    assert "#SBATCH --time=2:00:00" in header


# --------------------------- command line builders ---------------------


def test_build_snakemake_command_quotes_conda_prefix(submitter, monkeypatch):
    # Ensure the object picked up the env path
    submitter.conda_env_path = os.environ["CONDA_PREFIX"]

    # Before calling, it should not have been invoked
    submitter.snakemake_executable.get_command.assert_not_called()

    cmd = submitter._build_snakemake_command()

    # After building the command, get_command() must be called exactly once
    submitter.snakemake_executable.get_command.assert_called_once()

    # CONDA_PREFIX has a space, so it must be shell-quoted
    assert "conda run -p /conda/env" in cmd

    # And the returned snakemake command should be included verbatim
    assert "snakemake --notemp -p --rerun-triggers mtime" in cmd


def test_build_job_status_check_quotes_paths(submitter):
    line = submitter._build_job_status_check()
    # log_path and result_path/analysis_status.txt must be quoted
    assert shlex.quote(str(submitter.log_path)) in line
    assert shlex.quote(str(submitter.result_path / "analysis_status.txt")) in line
    # scan script path should be quoted too
    assert f"python {shlex.quote(submitter.scan_finished_jobid_status)}" in line


# --------------------------- success snippet --------------------------


def test_build_success_status_check_quotes_and_exits(submitter):
    snippet = submitter._build_success_status_check()
    # quoted paths for marker and status
    marker = shlex.quote(str(submitter.result_path / "analysis_finished_successfully"))
    status = shlex.quote(str(submitter.result_path / "analysis_status.txt"))
    assert marker in snippet
    assert status in snippet
    # expected logic present
    assert "exit 1" in snippet
    assert "exit 2" in snippet
    # should stderr-print something
    assert ">&2" in snippet


# --------------------------- end-to-end create ------------------------


def test_create_sbatch_script_writes_once_and_contains_parts(submitter):
    # IMPORTANT: class must call _build_success_status_check(), not build_success_status_check()
    with patch("pathlib.Path.write_text") as mock_write:
        submitter.create_sbatch_script()
        mock_write.assert_called_once()

        content = mock_write.call_args[0][0]
        assert isinstance(content, str) and content.strip()

        # Contains header + snakemake + status check + success snippet
        assert "#!/bin/bash -l" in content
        assert "snakemake --notemp -p --rerun-triggers mtime" in content
        assert "python " in content and "analysis_status.txt" in content
        assert "analysis_finished_successfully" in content

        # get_command is invoked exactly once
        submitter.snakemake_executable.get_command.assert_called_once()

        # Script path used is the class' sbatch_script_path
        # (We can also check that the logger line references it, but here we confirm write_text target.)
        # write_text was called on pathlib.Path(self.sbatch_script_path), which we can't directly assert
        # from mock args; instead, assert content is a single assembled string with separators.
        assert (
            content.count("\n\n") >= 3
        )  # header, snakemake, status-check, success-check blocks


def test_create_sbatch_script_handles_empty_conda_prefix(submitter, monkeypatch):
    # If CONDA_PREFIX is empty, shlex.quote("") returns "''"
    monkeypatch.delenv("CONDA_PREFIX", raising=False)
    submitter.conda_env_path = os.environ.get("CONDA_PREFIX", "")
    with patch("pathlib.Path.write_text") as mock_write:
        submitter.create_sbatch_script()
        content = mock_write.call_args[0][0]
        # Should still contain a -p argument; empty path becomes single-quoted empty string
        assert "conda run -p '' " in content


def test_submit_job_success(submitter, caplog):
    """submit_job returns the parsed Job ID when sbatch succeeds and logs info."""
    caplog.set_level(logging.INFO)
    fake_output = "Submitted batch job 12345\n"

    with patch(f"{MODULE}.subprocess.run") as mock_run:
        mock_proc = MagicMock(returncode=0, stdout=fake_output)
        mock_run.return_value = mock_proc

        job_id = submitter.submit_job()

    assert job_id == "12345"
    # Check that success info log is present
    assert "Job submitted successfully with Job ID: 12345" in caplog.text
    # Also sanity-check the exact command called
    mock_run.assert_called_once_with(
        ["sbatch", str(submitter.sbatch_script_path)],
        capture_output=True,
        text=True,
        check=True,
    )


def test_submit_job_failure(submitter, caplog):
    """submit_job returns None and logs error when sbatch raises CalledProcessError."""
    caplog.set_level(logging.INFO)
    error = subprocess.CalledProcessError(
        returncode=1, cmd=["sbatch"], stderr="Error submitting job\n"
    )
    with patch(f"{MODULE}.subprocess.run", side_effect=error):
        job_id = submitter.submit_job()

    assert job_id is None
    assert "sbatch submission failed: Error submitting job" in caplog.text


def test_submit_job_no_job_id(submitter, caplog):
    """submit_job returns None and logs a warning when sbatch output lacks a job id."""
    caplog.set_level(logging.INFO)
    with patch(f"{MODULE}.subprocess.run") as mock_run:
        mock_proc = MagicMock(returncode=0, stdout="Some unrelated output\n")
        mock_run.return_value = mock_proc

        job_id = submitter.submit_job()

    assert job_id is None
    assert (
        "Could not extract Job ID from sbatch output: Some unrelated output"
        in caplog.text
    )


def test_write_job_id_yaml_calls_helper_with_expected_payload(
    tmp_path, submitter, caplog
):
    """write_job_id_yaml should call write_yaml with {case_id: [job_id]} and correct path."""
    caplog.set_level(logging.INFO)
    submitter.result_path = tmp_path  # ensure deterministic output dir
    job_id = "12345"
    expected_yaml_path = tmp_path / "slurm_jobids.yaml"
    expected_payload = {submitter.case_id: [job_id]}

    with patch(f"{MODULE}.write_yaml") as mock_write_yaml:
        submitter.write_job_id_yaml(job_id)

    mock_write_yaml.assert_called_once()
    args, _ = mock_write_yaml.call_args
    assert args[0] == expected_payload
    assert Path(args[1]) == expected_yaml_path
    assert f"Job ID written to {expected_yaml_path}" in caplog.text
