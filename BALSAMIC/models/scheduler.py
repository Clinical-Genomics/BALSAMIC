"""Scheduler models."""
import logging
import subprocess
import sys
from pathlib import Path
from re import Match, search
from typing import Any, Dict, List, Optional

from pydantic import BaseModel, DirectoryPath, Field, FilePath, field_validator
from pydantic_core.core_schema import ValidationInfo

from BALSAMIC.constants.cluster import QOS, ClusterMailType, ClusterProfile
from BALSAMIC.utils.utils import remove_unnecessary_spaces

LOG = logging.getLogger(__name__)


class Scheduler(BaseModel):
    """
    Scheduler model handling cluster job submissions.

    Attributes:
        account (str)                             : Cluster account to run jobs.
        benchmark (Optional[bool])                : Flag to profile slurm jobs.
        case_id (str)                             : Case identifier.
        dependencies (Optional[List[str]])        : List of job dependencies.
        job_properties (Dict[str, Any])           : Job properties defined in a snakemake jobscript.
        job_script (FilePath)                     : Snakemake job script path.
        log_dir (DirectoryPath)                   : Logging directory.
        mail_type (Optional[ClusterMailType])     : Email type triggering job status notifications.
        mail_user (Optional[str])                 : User email to receive job status notifications.
        profile (ClusterProfile)                  : Cluster profile to submit jobs.
        profiling_interval (Optional[int])        : Sampling interval for a profiling type.
        profiling_type (Optional[str])            : Collected data types.
        qos (Optional[QOS])                       : QOS for sbatch jobs.
    """

    account: str
    benchmark: Optional[bool] = False
    case_id: str
    dependencies: Optional[List[str]] = None
    job_properties: Dict[str, Any]
    job_script: FilePath
    log_dir: DirectoryPath
    mail_type: Optional[ClusterMailType] = Field(default=None, validate_default=True)
    mail_user: Optional[str] = Field(default=None, validate_default=True)
    profile: ClusterProfile
    profiling_interval: int = 10
    profiling_type: str = "task"
    qos: Optional[QOS] = QOS.LOW

    @field_validator("account")
    def get_account_option(cls, account: str) -> str:
        """Return string representation of the account option."""
        return f"--account {account}"

    @field_validator("mail_type")
    def get_mail_type_option(
        cls, mail_type: Optional[ClusterMailType], info: ValidationInfo
    ) -> str:
        """Return string representation of the mail_type option."""
        if mail_type:
            return f"--mail-type {mail_type}"
        cluster_mail_type: Optional[ClusterMailType] = (
            info.data["job_properties"]["cluster"].get("mail_type")
            if info.data.get("job_properties")
            else None
        )
        if cluster_mail_type:
            return f"--mail-type {cluster_mail_type}"
        return ""

    @field_validator("mail_user")
    def get_mail_user_option(cls, mail_user: Optional[str]) -> str:
        """Return string representation of the mail_user option."""
        if mail_user:
            return f"--mail-user {mail_user}"
        return ""

    @field_validator("qos")
    def get_qos_option(cls, qos: Optional[QOS]) -> str:
        """Return string representation of the mail_user option."""
        if qos:
            return f"--qos {qos}"
        return ""

    def get_dependency_option(self) -> str:
        """Return string representation of the dependency option."""
        if self.dependencies:
            dependencies: str = ",".join(
                [job for job in self.dependencies if job.isdigit()]
            )
            return f"--dependency afterok:{dependencies}"
        return ""

    def get_error_option(self) -> str:
        """Return the standard error file path."""
        return f"--error {Path(self.log_dir, f'{Path(self.job_script).name}.%j.err').as_posix()}"

    def get_output_option(self) -> str:
        """Return the standard output file path."""
        return f"--output {Path(self.log_dir, f'{Path(self.job_script).name}.%j.out').as_posix()}"

    def get_profile_option(self) -> str:
        """Return string representation of the slurm profile option."""
        if self.benchmark and self.profile == ClusterProfile.SLURM:
            return f"--profile {self.profiling_type}"
        return ""

    def get_acctg_freq_option(self) -> str:
        """Return string representation of the profiling sampling intervals in seconds option."""
        if self.benchmark and self.profile == ClusterProfile.SLURM:
            return f"--acctg-freq {self.profiling_type}={self.profiling_interval}"
        return ""

    def get_ntasks_option(self) -> str:
        """Return the maximum number of tasks for allocation."""
        ntasks: str = self.job_properties["cluster"].get("n")
        if ntasks:
            return f"--ntasks {ntasks}"
        return ""

    def get_memory_option(self) -> str:
        """Return the maximum memory allocation."""
        mem: str = self.job_properties["cluster"].get("mem")
        if mem:
            return f"--mem {mem}MB"
        return ""

    def get_time_option(self) -> str:
        """Return the allocation time."""
        time: str = self.job_properties["cluster"].get("time")
        if time:
            return f"--time {time}"
        return ""

    def get_partition_option(self) -> str:
        """Return the specific partition for the resource allocation."""
        partition: str = self.job_properties["cluster"].get("partition")
        if partition:
            return f"--partition {partition}"
        return ""

    @staticmethod
    def get_job_id_from_stdout(stdout: str) -> str:
        """Return job ID from the standard output."""
        job_id_match: Match[str] = search("Submitted batch job (\d+)", stdout)
        if job_id_match:
            job_id: str = job_id_match.group(1)
            LOG.info(f"Submitted job with ID: {job_id}")
            return job_id
        LOG.error("Failed to extract job ID from the submission result")
        raise ValueError

    def write_job_log_data(self, job_id: str, command: str) -> None:
        """Write accounting information for jobs."""
        log_path: Path = Path(self.log_dir, f"{self.case_id}.sacct")
        with open(log_path, "a") as file:
            file.write(f"{job_id},{command}\n")

    def get_command(self) -> str:
        """Return the command to submit a specific job to the cluster."""
        command: str = (
            f"sbatch "
            f"{self.account} "
            f"{self.mail_type} "
            f"{self.mail_user} "
            f"{self.qos} "
            f"{self.get_dependency_option()} "
            f"{self.get_error_option()} "
            f"{self.get_output_option()} "
            f"{self.get_profile_option()} "
            f"{self.get_acctg_freq_option()} "
            f"{self.get_ntasks_option()} "
            f"{self.get_time_option()} "
            f"{self.get_memory_option()} "
            f"{self.get_partition_option()} "
            f"{self.job_script} "
        )
        return remove_unnecessary_spaces(command)

    def submit_job(self) -> None:
        """Submit a job to the cluster."""
        cluster_command: str = self.get_command()
        try:
            result: subprocess.CompletedProcess = subprocess.run(
                cluster_command,
                check=True,
                shell=True,
                stdout=subprocess.PIPE,
                text=True,
            )
            job_id: str = self.get_job_id_from_stdout(result.stdout)
            self.write_job_log_data(job_id=job_id, command=cluster_command)
            # Send job ID to stdout for dependency parsing
            print(job_id, file=sys.stdout)
        except Exception:
            LOG.error(f"Failed to submit: {cluster_command}")
            raise
