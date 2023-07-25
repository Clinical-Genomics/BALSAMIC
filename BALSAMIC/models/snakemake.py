"""Snakemake related models."""
import sys
from pathlib import Path
from typing import Dict, Optional, List

from pydantic import BaseModel, FilePath, DirectoryPath, validator

from BALSAMIC.constants.analysis import RunMode
from BALSAMIC.constants.cluster import ClusterMailType, QOS, ClusterProfile, MAX_JOBS
from BALSAMIC.constants.paths import SCHEDULER_PATH
from BALSAMIC.utils.utils import remove_unnecessary_spaces


class SingularityBindPath(BaseModel):
    """Singularity binding path model.

    Attributes:
        source (Path)      : Path to the file or directory on the host system.
        destination (Path) : Path inside the container where the source will be mounted.
    """

    source: Path
    destination: Path


class Snakemake(BaseModel):
    """Snakemake command building model.

    Attributes:
        account (Optional[str])                           : Scheduler account.
        benchmark (Optional[bool])                        : Slurm jobs profiling option.
        case_id (str)                                     : Analysis case name.
        cluster_config_path (Optional[FilePath])          : Cluster configuration file path.
        config_path (FilePath)                            : Sample configuration file.
        disable_variant_caller (Optional[str])            : Disable variant caller.
        dragen (Optional[bool])                           : FLag for enabling or disabling Dragen suite.
        force (bool)                                      : Force snakemake execution.
        log_dir (Optional[DirectoryPath])                 : Logging directory.
        mail_type (Optional[ClusterMailType])             : Email type triggering job status notifications.
        mail_user (Optional[str])                         : User email to receive job status notifications.
        profile (Optional[ClusterProfile])                : Cluster profile to submit jobs.
        qos (Optional[QOS])                               : QOS for sbatch jobs.
        quiet (Optional[bool])                            : Quiet mode for snakemake.
        report_path (Optional[Path])                      : Snakemake generated report path.
        result_dir (Optional[DirectoryPath])              : Analysis output directory.
        run_analysis (bool)                               : Flag to run the actual analysis.
        run_mode (RunMode)                                : Cluster run mode to execute analysis.
        script_dir (Optional[DirectoryPath])              : Cluster profile scripts directory.
        singularity_bind_paths (Optional[Dict[str, str]]) : Singularity source and destination bind paths.
        snakefile (FilePath)                              : Snakemake rule configuration file.
        snakemake_options (Optional[List[str]])           : Snakemake command additional options.
        working_dir (Path)                                : Snakemake working directory.
    """

    account: Optional[str]
    benchmark: bool = False
    case_id: str
    cluster_config_path: Optional[FilePath]
    config_path: FilePath
    disable_variant_caller: Optional[str]
    dragen: bool = False
    force: bool = False
    log_dir: Optional[DirectoryPath]
    mail_type: Optional[ClusterMailType]
    mail_user: Optional[str]
    profile: Optional[ClusterProfile]
    qos: Optional[QOS]
    quiet: bool = False
    report_path: Optional[Path]
    result_dir: Optional[DirectoryPath]
    run_analysis: bool = False
    run_mode: RunMode
    script_dir: Optional[DirectoryPath]
    singularity_bind_paths: Optional[List[SingularityBindPath]]
    snakefile: FilePath
    snakemake_options: Optional[List[str]]
    working_dir: Path

    @validator("disable_variant_caller", always=True)
    def get_disable_variant_caller_option(cls, disable_variant_caller: str) -> str:
        """Return string representation of the disable_variant_caller option."""
        if disable_variant_caller:
            return f"disable_variant_caller={disable_variant_caller}"
        return ""

    @validator("mail_user", always=True)
    def get_mail_user_option(cls, mail_user: Optional[str]) -> str:
        """Return string representation of the mail_user option."""
        if mail_user:
            return f"--mail-user {mail_user}"
        return ""

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

    def get_mail_type_option(self) -> str:
        """Return string representation of the mail_type option."""
        if self.mail_type:
            return f"--mail-type {self.mail_type.value}"
        return ""

    def get_quiet_flag(self) -> str:
        """Return string representation of the quiet flag."""
        if self.quiet:
            return "--quiet"
        return ""

    def get_report_path_option(self) -> str:
        """Return string representation of the report_path option."""
        if self.report_path:
            return f"--report {self.report_path.as_posix()}"
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
            return f"--use-singularity --singularity-args ' --cleanenv {' '.join(bind_options)}'"
        return ""

    def get_slurm_profiler_option(self) -> str:
        """Return string representation of the slurm profiler option."""
        if self.benchmark and self.profile == ClusterProfile.SLURM:
            return "--slurm-profiler task"
        return ""

    def get_snakemake_options_command(self) -> str:
        """Return string representation of the additional Snakemake options."""
        if self.snakemake_options:
            return " ".join(self.snakemake_options)
        return ""

    def get_snakemake_command(self) -> str:
        """Return Snakemake command to be submitted."""
        snakemake_command: str = (
            f"snakemake --notemp -p "
            f"--directory {self.working_dir.as_posix()} "
            f"--snakefile {self.snakefile.as_posix()} "
            f"--configfiles {self.config_path.as_posix()} "
            f"{self.cluster_config_path.as_posix()} "
            f"{self.get_singularity_bind_paths_option()} "
            f"{self.get_quiet_flag()} "
            f"{self.get_force_flag()} "
            f"{self.get_run_analysis_flag()} "
            f"{self.get_snakemake_cluster_options()} "
            f"{self.get_report_path_option()} "
            f"{self.get_snakemake_config_options()} "
            f"{self.get_snakemake_options_command()}"
        )
        return remove_unnecessary_spaces(snakemake_command)

    def get_snakemake_config_options(self) -> str:
        """Return Snakemake config options to be submitted."""
        return remove_unnecessary_spaces(
            f"--config {self.disable_variant_caller} {self.get_dragen_flag()}"
        )

    def get_snakemake_cluster_options(self) -> str:
        """Return Snakemake cluster options to be submitted."""
        if self.run_mode == RunMode.CLUSTER:
            snakemake_cluster_options: str = (
                f"--immediate-submit -j {MAX_JOBS} "
                f"--jobname BALSAMIC.{self.case_id}.{{rulename}}.{{jobid}}.sh "
                f"--cluster-config {self.cluster_config_path.as_posix()} "
                f"--cluster {self.get_cluster_submit_command()}"
            )
            return remove_unnecessary_spaces(snakemake_cluster_options)
        return ""

    def get_cluster_submit_command(self) -> str:
        """Get cluster command to be submitted by Snakemake."""
        cluster_submit_command: str = (
            f"'{sys.executable} {SCHEDULER_PATH.as_posix()} "
            f"--sample-config {self.config_path.as_posix()} "
            f"--profile {self.profile.value} "
            f"--account {self.account} "
            f"--qos {self.qos.value} "
            f"--log-dir {self.log_dir.as_posix()} "
            f"--script-dir {self.script_dir.as_posix()} "
            f"--result-dir {self.result_dir.as_posix()} "
            f"{self.get_slurm_profiler_option()} "
            f"{self.mail_user} "
            f"{self.get_mail_type_option()} "
            "{dependencies} '"
        )
        return remove_unnecessary_spaces(cluster_submit_command)
