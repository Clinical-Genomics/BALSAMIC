import os
import re
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
        self.max_run_hours = max_run_hours
        self.snakemake_executable = snakemake_executable
        self.log = logger

        self.conda_env_path = os.environ.get("CONDA_PREFIX", "")
        self.sbatch_script_path = self.script_path / "BALSAMIC_snakemake_submit.sh"

    def create_sbatch_script(self) -> None:
        self.log.info("Creating sbatch script to submit jobs.")
        self.log.info(f"Using conda environment: {self.conda_env_path}")

        sbatch_header = textwrap.dedent(
            f"""\
            #!/bin/bash -l
            #SBATCH --account={self.account}
            #SBATCH --job-name=BALSAMIC_snakemake_submit.{self.case_id}.%j
            #SBATCH --output={self.log_path}/BALSAMIC_snakemake_submit.{self.case_id}.%j.out
            #SBATCH --error={self.log_path}/BALSAMIC_snakemake_submit.{self.case_id}.%j.err
            #SBATCH --ntasks=1
            #SBATCH --mem=5G
            #SBATCH --time={self.max_run_hours}:00:00
            #SBATCH --qos={self.qos}
            #SBATCH --cpus-per-task=1
        """
        )

        # Run snakemake workflow
        sbatch_command = f"\nconda run -p {self.conda_env_path} {self.snakemake_executable.get_command()}\n"

        # Check the status of submitted jobs
        job_status_check = f"\nconda run -p {self.conda_env_path} python {self.scan_finished_jobid_status} {self.log_path} --output {self.result_path}/analysis_status.txt\n"

        # Check the final success status of the workflow
        success_status_check = textwrap.dedent(
            f"""\n
            if [[ -f "{self.result_path}/analysis_status.txt" ]]; then
                STATUS=$(cat "{self.result_path}/analysis_status.txt")
                echo "Snakemake analysis status: $STATUS"
                if [[ "$STATUS" != "SUCCESS" ]]; then
                    echo "Analysis failed: $STATUS"
                    exit 1
                fi
            else
                echo "No status file found; assuming failure"
                exit 2
            fi \n
        """
        )

        full_script = (
            sbatch_header + sbatch_command + job_status_check + success_status_check
        )

        with open(self.sbatch_script_path, "w") as f:
            f.write(full_script)

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
