"""Snakemake related models."""
import sys
from pathlib import Path
from typing import List, Optional

from pydantic import BaseModel, DirectoryPath, FilePath

from BALSAMIC.constants.analysis import RunMode
from BALSAMIC.constants.cluster import MAX_JOBS, QOS
from BALSAMIC.utils.utils import remove_unnecessary_spaces


class SingularityBindPath(BaseModel):
    """Singularity binding path model.

    Attributes:
        source (Path)      : Path to the file or directory on the host system.
        destination (Path) : Path inside the container where the source will be mounted.
    """

    source: Path
    destination: Path


class SnakemakeExecutable(BaseModel):
    """Snakemake command building model.

    Attributes:
        account (Optional[str])                                      : Scheduler account.
        case_id (str)                                                : Analysis case name.
        config_path (FilePath)                                       : Sample configuration file.
        dragen (Optional[bool])                                      : Flag for enabling or disabling Dragen suite.
        force (bool)                                                 : Force snakemake execution.
        log_dir (Optional[DirectoryPath])                            : Logging directory.
        cluster_profile: Path                                        : Directory containing snakemake cluster profile
        cluster_job_status_script (FilePath)                         : Path to script for snakemake to parse more slurm job-statuses
        workflow_profile: Path                                       : Directory contianing snakemake workflow profile specifying rule resources
        qos (Optional[QOS])                                          : QOS for sbatch jobs.
        quiet (Optional[bool])                                       : Quiet mode for snakemake.
        run_analysis (bool)                                          : Flag to run the actual analysis.
        run_mode (RunMode)                                           : Cluster run mode to execute analysis.
        script_dir (Optional[DirectoryPath])                         : Cluster profile scripts directory.
        singularity_bind_paths (Optional[List[SingularityBindPath]]) : Singularity source and destination bind paths.
        snakefile (FilePath)                                         : Snakemake rule configuration file.
        snakemake_options (Optional[List[str]])                      : Snakemake command additional options.
        working_dir (Path)                                           : Snakemake working directory.

    """

    account: Optional[str] = None
    case_id: str
    config_path: FilePath
    dragen: bool = False
    force: bool = False
    log_dir: Optional[DirectoryPath] = None
    cluster_profile: Path
    cluster_job_status_script: FilePath
    workflow_profile: Path
    qos: Optional[QOS] = None
    quiet: bool = False
    run_analysis: bool = False
    run_mode: RunMode
    script_dir: Optional[DirectoryPath] = None
    singularity_bind_paths: Optional[List[SingularityBindPath]] = None
    snakefile: FilePath
    snakemake_options: Optional[List[str]] = None
    working_dir: Path

    def get_config_file_option(self) -> str:
        """Return string representation of the config file."""
        return f"--configfile {self.config_path.as_posix()}"

    def get_config_options(self) -> str:
        """Return Snakemake config options to be submitted."""
        return remove_unnecessary_spaces(
            f"{f'--config {self.get_dragen_flag()}' if self.get_dragen_flag() else ''} "
        )

    def get_cluster_status_script(self) -> str:
        """Return cluster-status argument."""
        return f'--cluster-status "python {self.cluster_job_status_script.as_posix()}"'

    def get_dragen_flag(self) -> str:
        """Return string representation of the dragen flag."""
        if self.dragen:
            return "dragen=True"
        return ""

    def get_force_flag(self) -> str:
        """Return string representation of the force flag."""
        if self.force:
            return "--forceall"
        return ""

    def get_quiet_flag(self) -> str:
        """Return string representation of the quiet flag."""
        if self.quiet:
            return "--quiet"
        return ""

    def get_run_analysis_flag(self) -> str:
        """Return string representation of the run_analysis flag."""
        if not self.run_analysis:
            return "--dryrun"
        return ""

    def get_singularity_bind_paths_option(self) -> str:
        """Return string representation of the singularity_bind_paths option."""
        if self.singularity_bind_paths:
            bind_options: List[str] = []
            for singularity_bind_path in self.singularity_bind_paths:
                bind_options.append(
                    f"--bind {singularity_bind_path.source.as_posix()}:{singularity_bind_path.destination.as_posix()}"
                )
            return f"--use-singularity --singularity-args '--cleanenv {' '.join(bind_options)}'"
        return ""

    def get_snakemake_options_command(self) -> str:
        """Return string representation of the additional Snakemake options."""
        if self.snakemake_options:
            return " ".join(self.snakemake_options)
        return ""

    def get_command(self) -> str:
        """Return Snakemake command to be submitted."""
        snakemake_command: str = (
            f"snakemake --notemp -p --rerun-triggers mtime "
            f"--directory {self.working_dir.as_posix()} "
            f"--snakefile {self.snakefile.as_posix()} "
            f"{self.get_config_file_option()} "
            f"{self.get_config_options()} "
            f"{self.get_singularity_bind_paths_option()} "
            f"{self.get_quiet_flag()} "
            f"{self.get_force_flag()} "
            f"--slurm-logdir {self.log_dir} "
            f"{self.get_run_analysis_flag()} "
            f"{self.get_snakemake_cluster_options()} "
            "--executor slurm "
            f"{self.get_snakemake_options_command()}"
        )
        return remove_unnecessary_spaces(snakemake_command)

    def get_snakemake_cluster_options(self) -> str:
        """Return Snakemake cluster options to be submitted."""
        if self.run_mode == RunMode.CLUSTER:
            snakemake_cluster_options: str = (
                f"-j {MAX_JOBS} "
                f"--profile {self.workflow_profile} "
                f"--jobname BALSAMIC.{self.case_id}.{{rulename}}.{{jobid}} "
                f"--default-resources slurm_extra=\"--qos={self.qos} --output {self.log_dir}/%x.%j.out --error {self.log_dir}/%x.%j.err\" "
                f"slurm_partition=core slurm_account={self.account} "
                f"--slurm-keep-successful-logs"
            )
            return remove_unnecessary_spaces(snakemake_cluster_options)
        return "--default-resources mem_mb=32000 threads=8"
