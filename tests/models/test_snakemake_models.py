"""Tests for the Snakemake model related methods."""
import copy
import sys
from pathlib import Path
from typing import Any, Dict

import pytest
from BALSAMIC.constants.cluster import MAX_JOBS, QOS, ClusterAccount, ClusterProfile
from BALSAMIC.constants.paths import IMMEDIATE_SUBMIT_PATH
from BALSAMIC.models.snakemake import SingularityBindPath, SnakemakeExecutable
from pydantic import ValidationError


def test_singularity_bind_path_model(singularity_bind_path_data: Dict[str, Path]):
    """Test singularity bind path model initialisation."""

    # GIVEN singularity bind path model data

    # WHEN initialising the model
    singularity_bind_path_model: SingularityBindPath = SingularityBindPath(
        **singularity_bind_path_data
    )

    # THEN the model should have been correctly built
    assert dict(singularity_bind_path_model) == singularity_bind_path_data


def test_snakemake_model(
    snakemake_executable_data: Dict[str, Any],
    snakemake_executable_validated_data: Dict[str, Any],
):
    """Test snakemake model initialisation."""

    # GIVEN a cluster ready snakemake data

    # WHEN initialising the model
    snakemake_model: SnakemakeExecutable = SnakemakeExecutable(
        **snakemake_executable_data
    )

    # THEN the model should have been correctly built
    assert dict(snakemake_model) == snakemake_executable_validated_data


def test_snakemake_model_empty():
    """Test Snakemake empty model initialisation."""

    # GIVEN no input for the snakemake model

    # WHEN initialising the model
    with pytest.raises(ValidationError):
        # THEN an empty model should raise a ValidationError
        SnakemakeExecutable()


def test_get_config_files_option(
    snakemake_executable: SnakemakeExecutable, reference_file: Path
):
    """Test formatting of the configuration files."""

    # GIVEN a snakemake executable model with a mocked config file

    # WHEN calling the method
    config_files_option: str = snakemake_executable.get_config_files_option()

    # THEN the expected format should be returned
    assert (
        config_files_option
        == f"--configfiles {reference_file.as_posix()} {reference_file.as_posix()}"
    )


def test_get_config_options(snakemake_executable: SnakemakeExecutable):
    """Test formatting of the snakemake config options."""

    # GIVEN a snakemake executable model disabling some variant callers

    # WHEN calling the method
    snakemake_config_options: str = snakemake_executable.get_config_options()

    # THEN the expected format should be returned
    assert snakemake_config_options == "--config disable_variant_caller=tnscope"


def test_get_dragen_flag(snakemake_executable: SnakemakeExecutable):
    """Test formatting of the dragen flag."""

    # GIVEN a snakemake executable model with a dragen flag
    snakemake_model: SnakemakeExecutable = copy.deepcopy(snakemake_executable)
    snakemake_model.dragen = True

    # WHEN calling the method
    dragen_flag: str = snakemake_model.get_dragen_flag()

    # THEN the expected format should be returned
    assert dragen_flag == "dragen=True"


def test_get_force_flag(snakemake_executable: SnakemakeExecutable):
    """Test formatting of the force flag."""

    # GIVEN a snakemake executable model with a force flag
    snakemake_model: SnakemakeExecutable = copy.deepcopy(snakemake_executable)
    snakemake_model.force = True

    # WHEN calling the method
    force_flag: str = snakemake_model.get_force_flag()

    # THEN the expected format should be returned
    assert force_flag == "--forceall"


def test_quiet_flag(snakemake_executable: SnakemakeExecutable):
    """Test formatting of the quiet flag."""

    # GIVEN a snakemake executable model with a quiet option

    # WHEN calling the method
    quiet_flag: str = snakemake_executable.get_quiet_flag()

    # THEN the expected format should be returned
    assert quiet_flag == "--quiet"


def test_get_run_analysis_flag(snakemake_executable: SnakemakeExecutable):
    """Test formatting of the run analysis flag."""

    # GIVEN a snakemake executable model with a dry run option
    snakemake_model: SnakemakeExecutable = copy.deepcopy(snakemake_executable)
    snakemake_model.run_analysis = False

    # WHEN calling the method
    run_analysis_flag: str = snakemake_model.get_run_analysis_flag()

    # THEN the expected format should be returned
    assert run_analysis_flag == "--dryrun"


def test_get_singularity_bind_paths_option(
    reference_file: Path,
    session_tmp_path: Path,
    snakemake_executable: SnakemakeExecutable,
):
    """Test formatting of the singularity bind paths."""

    # GIVEN a snakemake executable model with multiple binding paths
    snakemake_model: SnakemakeExecutable = copy.deepcopy(snakemake_executable)
    snakemake_model.singularity_bind_paths.append(
        SingularityBindPath(source=reference_file, destination=reference_file)
    )

    # WHEN calling the method
    singularity_bind_paths_option: str = (
        snakemake_model.get_singularity_bind_paths_option()
    )

    # THEN the expected format should be returned
    assert (
        singularity_bind_paths_option
        == f"--use-singularity --singularity-args '--cleanenv --bind {session_tmp_path.as_posix()}:/ "
        f"--bind {reference_file.as_posix()}:{reference_file.as_posix()}'"
    )


def test_get_snakemake_options_command(snakemake_executable: SnakemakeExecutable):
    """Test formatting of the snakemake options command."""

    # GIVEN a snakemake executable model with additional snakemake options command

    # WHEN calling the method
    snakemake_options_command: str = (
        snakemake_executable.get_snakemake_options_command()
    )

    # THEN the expected format should be returned
    assert snakemake_options_command == "--cores 36"


def test_get_cluster_submit_command(
    case_id_tumor_only: str,
    mail_user_option: str,
    session_tmp_path: Path,
    snakemake_executable: SnakemakeExecutable,
):
    """Test formatting of the cluster submit command."""

    # GIVEN a snakemake executable model with working environment paths

    # WHEN calling the method
    snakemake_cluster_submit_command: str = (
        snakemake_executable.get_cluster_submit_command()
    )

    # THEN the expected format should be returned
    assert snakemake_cluster_submit_command == (
        f"'{sys.executable} {IMMEDIATE_SUBMIT_PATH.as_posix()} --account {ClusterAccount.DEVELOPMENT} "
        f"--log-dir {session_tmp_path} --mail-user {mail_user_option} --profile {ClusterProfile.SLURM} "
        f"--qos {QOS.HIGH} --script-dir {session_tmp_path} {case_id_tumor_only} {{dependencies}}'"
    )


def test_get_snakemake_cluster_options(
    case_id_tumor_only: str,
    mail_user_option: str,
    reference_file: Path,
    session_tmp_path: Path,
    snakemake_executable: SnakemakeExecutable,
):
    """Test formatting of the snakemake cluster options."""

    # GIVEN a snakemake executable model with working environment paths

    # WHEN calling the method
    snakemake_cluster_options: str = (
        snakemake_executable.get_snakemake_cluster_options()
    )

    # THEN the expected format should be returned
    assert (
        snakemake_cluster_options
        == f"--immediate-submit -j {MAX_JOBS} --jobname BALSAMIC.{case_id_tumor_only}.{{rulename}}.{{jobid}}.sh "
        f"--cluster-config {reference_file.as_posix()} --cluster '{sys.executable} {IMMEDIATE_SUBMIT_PATH.as_posix()} "
        f"--account {ClusterAccount.DEVELOPMENT} --log-dir {session_tmp_path} --mail-user {mail_user_option} "
        f"--profile {ClusterProfile.SLURM} --qos {QOS.HIGH} --script-dir {session_tmp_path} {case_id_tumor_only} "
        "{dependencies}'"
    )


def test_get_snakemake_command(
    case_id_tumor_only: str,
    mail_user_option: str,
    reference_file: Path,
    session_tmp_path: Path,
    snakemake_executable: SnakemakeExecutable,
):
    """Test retrieval of the snakemake command to be submitted to Slurm."""

    # GIVEN a snakemake executable model with working environment paths

    # WHEN calling the method
    snakemake_command: str = snakemake_executable.get_command()

    # THEN the expected format should be returned
    assert (
        snakemake_command
        == f"snakemake --notemp -p --rerun-trigger mtime --directory {session_tmp_path.as_posix()} "
        f"--snakefile {reference_file.as_posix()} "
        f"--configfiles {reference_file.as_posix()} {reference_file.as_posix()} "
        f"--use-singularity --singularity-args '--cleanenv --bind {session_tmp_path.as_posix()}:/' --quiet "
        f"--immediate-submit -j {MAX_JOBS} --jobname BALSAMIC.{case_id_tumor_only}.{{rulename}}.{{jobid}}.sh "
        f"--cluster-config {reference_file.as_posix()} --cluster '{sys.executable} {IMMEDIATE_SUBMIT_PATH.as_posix()} "
        f"--account {ClusterAccount.DEVELOPMENT} --log-dir {session_tmp_path} --mail-user {mail_user_option} "
        f"--profile {ClusterProfile.SLURM} --qos {QOS.HIGH} --script-dir {session_tmp_path} {case_id_tumor_only} "
        "{dependencies}' --config disable_variant_caller=tnscope --cores 36"
    )
