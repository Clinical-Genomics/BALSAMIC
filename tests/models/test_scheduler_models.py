"""Tests for the Scheduler model related methods."""
import logging
import subprocess
from pathlib import Path
from typing import Any, Dict
from unittest import mock

import pytest
from _pytest.logging import LogCaptureFixture
from pydantic import ValidationError

from BALSAMIC.constants.constants import EXIT_SUCCESS
from BALSAMIC.models.scheduler import Scheduler


def test_scheduler_model(
    scheduler_data: Dict[str, Any], scheduler_validated_data: Dict[str, Any]
):
    """Test immediate submit scheduler model initialisation."""

    # GIVEN a scheduler dictionary data

    # WHEN initialising the model
    scheduler_model: Scheduler = Scheduler(**scheduler_data)

    # THEN the model should have been correctly built
    assert scheduler_model.model_dump() == scheduler_validated_data


def test_scheduler_model_empty():
    """Test scheduler empty model initialisation."""

    # GIVEN no input for the scheduler model

    # WHEN initialising the model
    with pytest.raises(ValidationError):
        # THEN an empty model should raise a ValidationError
        Scheduler()


def test_get_error_option(scheduler_model: Scheduler):
    """Test get error option for cluster submission."""

    # GIVEN a scheduler model

    # WHEN retrieving the formatted option
    error_option: str = scheduler_model.get_error_option()

    # THEN the error option should be appended to the command
    assert "--error" in error_option


def test_get_output_option(scheduler_model: Scheduler):
    """Test parsing of the stdout option."""

    # GIVEN a scheduler model

    # WHEN getting the output option from the scheduler
    output_option: str = scheduler_model.get_output_option()

    # THEN the output option should be appended to the command
    assert "--output" in output_option


def test_get_profile_option(scheduler_model: Scheduler):
    """Test get scheduler profile option."""

    # GIVEN a scheduler model with a benchmarking option
    scheduler_model.benchmark = True

    # WHEN retrieving the profile option
    profile_option: str = scheduler_model.get_profile_option()

    # THEN the profile option should be appended to the command
    assert "--profile" in profile_option


def test_get_acctg_freq_option(scheduler_model: Scheduler):
    """Test get profiling sampling intervals."""

    # GIVEN a scheduler model with a benchmarking option
    scheduler_model.benchmark = True

    # WHEN getting the profiling option
    freq_option: str = scheduler_model.get_acctg_freq_option()

    # THEN the profiling option should be appended to the command
    assert "--acctg-freq" in freq_option


def test_get_ntasks_option(scheduler_model: Scheduler):
    """Test number of tasks for cluster resources allocation."""

    # GIVEN a scheduler model

    # WHEN getting the ntasks cluster resource option
    ntasks_option: str = scheduler_model.get_ntasks_option()

    # THEN the number of tasks option should be returned
    assert "--ntasks" in ntasks_option


def test_get_time_option(scheduler_model: Scheduler):
    """Test get time resource from job properties."""

    # GIVEN a scheduler model

    # WHEN getting the time cluster resource option
    time_option: str = scheduler_model.get_time_option()

    # THEN the number of tasks option should be returned
    assert "--time" in time_option


def test_get_partition_option(scheduler_model: Scheduler):
    """Test get partition from cluster properties."""

    # GIVEN a scheduler model

    # WHEN getting the partition cluster resource option
    partition_option: str = scheduler_model.get_partition_option()

    # THEN the partition option should be returned
    assert "--partition" in partition_option


def test_get_job_id_from_stdout(job_id: str, scheduler_model: Scheduler):
    """Test get job identifier from the scheduler standard output."""

    # GIVEN a scheduler model and a stdout with a job identifier
    stdout: str = f"Submitted batch job {job_id}"

    # WHEN retrieving the job id
    retrieved_job_id: str = scheduler_model.get_job_id_from_stdout(stdout)

    # THEN the expected job identifier should be returned
    assert retrieved_job_id == job_id


def test_get_job_id_from_stdout_error(
    job_id: str, scheduler_model: Scheduler, caplog: LogCaptureFixture
):
    """Test get job identifier from the scheduler incorrect standard output."""

    # GIVEN a scheduler model and an incorrect stdout with a job identifier
    stdout: str = "Submitted batch job no ID"

    # WHEN retrieving the job id
    with pytest.raises(ValueError):
        scheduler_model.get_job_id_from_stdout(stdout)

    # THEN the expected error should be catched
    assert "Failed to extract job ID from the submission result" in caplog.text


def test_write_job_log_data(job_id: str, scheduler_model: Scheduler):
    """Test writing job log data."""

    # GIVEN a scheduler model, a case id, a job id, and a scheduler command
    command: str = "sbatch --acount development --dependency 'afterok:00001'"

    # WHEN writing to the log file
    scheduler_model.write_job_log_data(job_id=job_id, command=command)

    # THEN the log file should have been created
    log_file: Path = Path(scheduler_model.log_dir, f"{scheduler_model.case_id}.sacct")
    assert log_file.is_file()
    with open(log_file, "r") as file:
        assert f"{job_id},{command}" in file.read()


def test_get_command(
    scheduler_model: Scheduler, scheduler_validated_data: Dict[str, Any]
):
    """Test scheduler command build."""

    # GIVEN a scheduler model and a validated data dictionary

    # WHEN getting the command to be submitted
    command: str = scheduler_model.get_command()

    # THEN the mandatory options should be appended to the command
    for option, value in scheduler_validated_data.items():
        if option not in ["benchmark", "case_id", "job_properties", "profile"]:
            assert str(value) in command


def test_submit_job(job_id: str, scheduler_model: Scheduler, caplog: LogCaptureFixture):
    """Test job submission to the cluster."""
    caplog.set_level(logging.INFO)

    # GIVEN a scheduler model and the expected standard output command
    stdout: str = f"Submitted batch job {job_id}"

    # WHEN submitting a cluster job
    with mock.patch.object(
        subprocess,
        "run",
        return_value=subprocess.CompletedProcess(
            args="", returncode=EXIT_SUCCESS, stdout=stdout, stderr=""
        ),
    ):
        scheduler_model.submit_job()

    # THEN the cluster job should be submitted and the log file created
    assert f"Submitted job with ID: {job_id}" in caplog.text
    assert Path(scheduler_model.log_dir, f"{scheduler_model.case_id}.sacct").is_file()


def test_submit_job_error(
    job_id: str, scheduler_model: Scheduler, caplog: LogCaptureFixture
):
    """Test job submission to the cluster when an exception is raised."""

    # GIVEN a scheduler model and the expected standard output command
    stdout: str = f"Submitted batch job {job_id}"

    # WHEN submitting a cluster job that has not been mocked
    with pytest.raises(Exception):
        scheduler_model.submit_job()

    # THEN the cluster job should fail to be submitted
    assert f"Failed to submit: {scheduler_model.get_command()}" in caplog.text
