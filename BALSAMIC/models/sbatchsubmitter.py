import os
import re
import textwrap
import subprocess
from pathlib import Path
from typing import Optional
from BALSAMIC.utils.io import write_json, write_yaml


class SbatchSubmitter:
    def __init__(
        self,
        case_id: str,
        script_path: Path,
        result_path: Path,
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

        sbatch_command = f"\nconda run -p {self.conda_env_path} {self.snakemake_executable.get_command()}\n"

        full_script = sbatch_header + sbatch_command

        with open(self.sbatch_script_path, "w") as f:
            f.write(full_script)

        self.log.info(f"Sbatch script written to: {self.sbatch_script_path}")

    def submit_job(self) -> Optional[str]:
        sbatch_command = f"sbatch {self.sbatch_script_path}"
        self.log.info(f"Submitting job with command: {sbatch_command}")

        result = subprocess.run(
            sbatch_command, shell=True, capture_output=True, text=True
        )

        if result.returncode == 0:
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
        else:
            self.log.error(f"sbatch submission failed: {result.stderr.strip()}")

        return None

    def write_job_id_yaml(self, job_id: str) -> None:
        yaml_path = self.result_path / "slurm_jobids.yaml"
        write_yaml({self.case_id: [job_id]}, yaml_path)
        self.log.info(f"Job ID written to {yaml_path}")
