"""Test Sbatch submitter class"""
from BALSAMIC.models.sbatchsubmitter import SbatchSubmitter
from unittest.mock import patch, MagicMock
import subprocess

MODULE = "BALSAMIC.models.sbatchsubmitter"


def test_create_sbatch_script(submitter):
    """create_sbatch_script should write exactly one script file via Path.write_text."""
    with patch("pathlib.Path.write_text") as mock_write_text:
        submitter.create_sbatch_script()
        mock_write_text.assert_called_once()
        # The script content should be a non-empty string
        (args, kwargs) = mock_write_text.call_args
        assert isinstance(args[0], str) and args[0].strip()


def test_submit_job_success(submitter):
    """submit_job returns the parsed Job ID when sbatch succeeds."""
    fake_output = "Submitted batch job 12345"
    with patch(f"{MODULE}.subprocess.run") as mock_run:
        mock_proc = MagicMock(returncode=0, stdout=fake_output)
        mock_run.return_value = mock_proc

        job_id = submitter.submit_job()

        assert job_id == "12345"
        submitter.log.info.assert_any_call(
            "Job submitted successfully with Job ID: 12345"
        )


def test_submit_job_failure(submitter):
    """submit_job returns None and logs error when sbatch raises CalledProcessError."""
    error = subprocess.CalledProcessError(
        returncode=1, cmd="sbatch", stderr="Error submitting job"
    )
    with patch(f"{MODULE}.subprocess.run", side_effect=error):
        job_id = submitter.submit_job()
        assert job_id is None
        submitter.log.error.assert_called_once()


def test_submit_job_no_job_id(submitter):
    """submit_job returns None and logs a warning when sbatch output lacks a job id."""
    with patch(f"{MODULE}.subprocess.run") as mock_run:
        mock_proc = MagicMock(returncode=0, stdout="Some unrelated output")
        mock_run.return_value = mock_proc

        job_id = submitter.submit_job()

        assert job_id is None
        submitter.log.warning.assert_called_once()


def test_write_job_id_yaml(tmp_path, submitter):
    """write_job_id_yaml writes a YAML file to result_path; ensure it lands on disk."""
    submitter.result_path = tmp_path  # ensure we write into a temp directory
    job_id = "12345"

    submitter.write_job_id_yaml(job_id)

    yaml_path = tmp_path / "slurm_jobids.yaml"
    assert yaml_path.is_file()
    # Optional: sanity-check that the job id appears in the file
    assert job_id in yaml_path.read_text()
