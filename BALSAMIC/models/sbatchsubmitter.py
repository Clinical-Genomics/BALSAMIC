import os
import re
import shlex
import textwrap
import subprocess
from pathlib import Path
from typing import Optional
from BALSAMIC.utils.io import write_yaml


class SbatchSubmitter:
    """SLURM job submission model for running a Snakemake workflow.

    Attributes:
        case_id (str)                  : Identifier for the analysis case.
        script_path (Path)            : Directory where the sbatch script will be created.
        result_path (Path)            : Directory where the job ID YAML file will be written.
        check_jobid_status_script (str): Python script for reporting failed or cancelled statuses of jobs in logdir
        log_path (Path)               : Directory where SLURM output and error logs will be written.
        account (str)                 : SLURM account to charge for the job.
        qos (str)                     : SLURM quality of service level.
        headjob_partition: Optional(str): Cluster partition for the headjob
        max_run_hours (int)          : Maximum allowed run time for the job, in hours.
        snakemake_executable         : Object representing the Snakemake command to be executed.
        logger                        : Logger instance for capturing logs.
        conda_env_path (str)         : Path to the active conda environment, from $CONDA_PREFIX.
        sbatch_script_path (Path)    : Path to the generated sbatch script file.
    """

    def __init__(
        self,
        case_id: str,
        script_path: Path,
        result_path: Path,
        scan_finished_jobid_status: str,
        log_path: Path,
        account: str,
        qos: str,
        headjob_partition: Optional[str],
        max_run_hours: int,
        snakemake_executable,
        logger,
    ):
        self.case_id = case_id
        self.script_path = script_path
        self.result_path = result_path
        self.scan_finished_jobid_status = scan_finished_jobid_status
        self.log_path = log_path
        self.account = account
        self.qos = qos
        self.headjob_partition = headjob_partition
        self.max_run_hours = max_run_hours
        self.snakemake_executable = snakemake_executable
        self.log = logger

        self.conda_env_path = os.environ.get("CONDA_PREFIX", "")
        self.sbatch_script_path = self.script_path / "BALSAMIC_snakemake_submit.sh"

    def _build_sbatch_header(self) -> str:
        """
        Construct the SBATCH header lines and return as a single string.

        Includes:
          - account, job-name, output/error paths, ntasks, mem, time, qos, cpus-per-task
          - optional partition if `self.headjob_partition` is set
        """
        lines = [
            "#!/bin/bash",
            f"#SBATCH --account={self.account}",
            f"#SBATCH --job-name=BALSAMIC_snakemake_submit.{self.case_id}.%j",
            f"#SBATCH --output={self.log_path}/BALSAMIC_snakemake_submit.{self.case_id}.%j.out",
            f"#SBATCH --error={self.log_path}/BALSAMIC_snakemake_submit.{self.case_id}.%j.err",
            "#SBATCH --ntasks=1",
            "#SBATCH --mem=5G",
            f"#SBATCH --time={self.max_run_hours}:00:00",
            f"#SBATCH --qos={self.qos}",
            "#SBATCH --cpus-per-task=1",
        ]
        if self.headjob_partition:
            lines.insert(2, f"#SBATCH --partition={self.headjob_partition}")
        return "\n".join(lines)

    def _build_snakemake_command(self) -> str:
        """
        Return the line that invokes Snakemake via conda.
        Uses shlex.quote for CONDA_PREFIX safety.
        """
        return (
            f"conda run -p {shlex.quote(self.conda_env_path)} "
            f"{self.snakemake_executable.get_command()}"
        )

    def _build_job_status_check(self) -> str:
        """
        Return the line that runs the job status checker via conda.
        Quotes all paths to be shell-safe.
        """
        return (
            f"conda run -p {shlex.quote(self.conda_env_path)} "
            f"python {shlex.quote(self.scan_finished_jobid_status)} "
            f"{shlex.quote(str(self.log_path))} --output "
            f"{shlex.quote(str(self.result_path / 'analysis_status.txt'))}"
        )

    def _build_success_status_check(
        self,
        success_marker: str = "analysis_finished_successfully",
        status_filename: str = "analysis_status.txt",
    ) -> str:
        """
        Build the bash snippet that verifies analysis success.

        Logic:
          - If <result_dir>/<success_marker> does NOT exist:
              - If <result_dir>/<status_filename> exists:
                  - Print a tag line and its contents to stderr, then exit 1
              - Else: print a generic error message and exit 2

        Returns:
            A dedented, stripped bash script fragment as a string.
        """
        # Quote paths for safe shell interpolation (handles spaces etc.)
        q_success = shlex.quote(str(Path(self.result_path, success_marker)))
        q_status = shlex.quote(str(Path(self.result_path, status_filename)))

        return textwrap.dedent(
            f"""
            if [[ ! -f {q_success} ]]; then
                if [[ -f {q_status} ]]; then
                    STATUS=$(cat {q_status})
                    echo "FROM ANALYSIS STATUS: {q_status}" >&2
                    echo "$STATUS" >&2
                    exit 1
                else
                    echo "No status file found; assuming error" >&2
                    exit 2
                fi
            fi
            """
        ).strip()

    def create_sbatch_script(self) -> None:
        """
        Generate and write an `sbatch` submission script for running a Snakemake workflow.
        """
        self.log.info("Creating sbatch script to submit jobs.")
        self.log.info(f"Using conda environment: {self.conda_env_path}")

        sbatch_header = self._build_sbatch_header()
        snakemake_cmd = self._build_snakemake_command()
        job_status_check = self._build_job_status_check()
        success_status_check = self._build_success_status_check()

        full_script = "\n\n".join(
            [sbatch_header, snakemake_cmd, job_status_check, success_status_check, ""]
        )

        Path(self.sbatch_script_path).write_text(full_script)
        self.log.info(f"Sbatch script written to: {self.sbatch_script_path}")

    def submit_job(self) -> Optional[str]:
        """Submit the generated sbatch script to the SLURM scheduler.

        Returns:
            Optional[str]: The SLURM job ID if the submission is successful, otherwise None.
        """
        command = ["sbatch", str(self.sbatch_script_path)]
        self.log.info(f"Submitting job with command: {' '.join(command)}")

        try:
            result = subprocess.run(
                command,
                capture_output=True,
                text=True,
                check=True,
            )
            output = result.stdout.strip()
            match = re.search(r"Submitted batch job (\d+)", output)
            if match:
                job_id = match.group(1)
                self.log.info(f"Job submitted successfully with Job ID: {job_id}")
                return job_id
            else:
                self.log.warning(
                    f"Could not extract Job ID from sbatch output: {output}"
                )
        except subprocess.CalledProcessError as e:
            self.log.error(f"sbatch submission failed: {e.stderr.strip()}")

        return None

    def write_job_id_yaml(self, job_id: str) -> None:
        """Write the submitted job ID to a YAML file.

        The file is saved at `case_id/analysis/slurm_jobids.yaml` and stores the job ID
        under the corresponding `case_id`.

        Args:
            job_id (str): The SLURM job ID to record.
        """
        yaml_path = self.result_path / "slurm_jobids.yaml"
        write_yaml({self.case_id: [job_id]}, yaml_path)
        self.log.info(f"Job ID written to {yaml_path}")
