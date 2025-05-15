"""Test Sbatch submitter class"""

from unittest.mock import patch, mock_open, MagicMock
from pathlib import Path


def test_create_sbatch_script(submitter):
    """Test that `create_sbatch_script` writes an sbatch script file.

    Ensures that the method attempts to open the target script file for writing
    and that content is written to it.
    """
    with patch("builtins.open", mock_open()) as mock_file:
        submitter.create_sbatch_script()
        # Check that the file was attempted to be written
        mock_file.assert_called_once_with(submitter.sbatch_script_path, "w")
        handle = mock_file()
        handle.write.assert_called()  # Basic check that something was written


def test_submit_job_success(submitter):
    """Test `submit_job` returns the job ID on successful sbatch submission.

    Simulates a valid SLURM response and verifies that the job ID is correctly parsed
    and logged.
    """
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
    """Test `submit_job` handles a failed sbatch call.

    Ensures that if sbatch returns a non-zero exit code, the method logs an error
    and returns None.
    """
    with patch("subprocess.run") as mock_run:
        mock_run.return_value.returncode = 1
        mock_run.return_value.stderr = "Error submitting job"
        job_id = submitter.submit_job()
        assert job_id is None
        submitter.log.error.assert_called_once()


def test_submit_job_no_job_id(submitter):
    """Test `submit_job` handles a successful sbatch call with unexpected output.

    Verifies that if sbatch output does not include a recognizable job ID,
    a warning is logged and None is returned.
    """
    with patch("subprocess.run") as mock_run:
        mock_run.return_value.returncode = 0
        mock_run.return_value.stdout = "Some unrelated output"
        job_id = submitter.submit_job()
        assert job_id is None
        submitter.log.warning.assert_called_once()


def test_write_job_id_yaml(submitter):
    """Test that `write_job_id_yaml` writes the job ID to a YAML file.

    Confirms that after calling the method, the expected YAML file exists
    in the result path.
    """
    job_id = "12345"
    submitter.write_job_id_yaml(job_id)
    assert Path(submitter.result_path, "slurm_jobids.yaml").is_file()
