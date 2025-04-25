"""Test Sbatch submitter class"""

from unittest.mock import patch, mock_open
from pathlib import Path


def test_create_sbatch_script(submitter):
    with patch("builtins.open", mock_open()) as mock_file:
        submitter.create_sbatch_script()
        # Check that the file was attempted to be written
        mock_file.assert_called_once_with(submitter.sbatch_script_path, "w")
        handle = mock_file()
        handle.write.assert_called()  # Basic check that something was written


def test_submit_job_success(submitter):
    fake_output = "Submitted batch job 12345"
    with patch("subprocess.run") as mock_run:
        mock_run.return_value.returncode = 0
        mock_run.return_value.stdout = fake_output
        job_id = submitter.submit_job()
        assert job_id == "12345"
        submitter.log.info.assert_any_call(
            "Job submitted successfully with Job ID: 12345"
        )


def test_submit_job_failure(submitter):
    with patch("subprocess.run") as mock_run:
        mock_run.return_value.returncode = 1
        mock_run.return_value.stderr = "Error submitting job"
        job_id = submitter.submit_job()
        assert job_id is None
        submitter.log.error.assert_called_once()


def test_submit_job_no_job_id(submitter):
    with patch("subprocess.run") as mock_run:
        mock_run.return_value.returncode = 0
        mock_run.return_value.stdout = "Some unrelated output"
        job_id = submitter.submit_job()
        assert job_id is None
        submitter.log.warning.assert_called_once()


def test_write_job_id_yaml(submitter):
    job_id = "12345"
    submitter.write_job_id_yaml(job_id)
    assert Path(submitter.result_path, "slurm_jobids.yaml").is_file()
