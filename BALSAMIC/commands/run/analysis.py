"""Balsamic run analysis CLI."""
import json
import logging
import os
import subprocess
import sys
from pathlib import Path
from typing import List

import click
import yaml

from BALSAMIC import __version__ as balsamic_version

from BALSAMIC.commands.options import (
    OPTION_CLUSTER_ACCOUNT,
    OPTION_CLUSTER_MAIL,
    OPTION_CLUSTER_MAIL_TYPE,
    OPTION_WORKFLOW_PROFILE,
    OPTION_CLUSTER_PROFILE,
    OPTION_CLUSTER_QOS,
    OPTION_CLUSTER_ENV,
    OPTION_DRAGEN,
    OPTION_FORCE_ALL,
    OPTION_QUIET,
    OPTION_RUN_ANALYSIS,
    OPTION_RUN_MODE,
    OPTION_RUN_INTERACTIVELY,
    OPTION_SAMPLE_CONFIG,
    OPTION_SNAKEFILE,
    OPTION_SNAKEMAKE_OPT,
)
from BALSAMIC.constants.analysis import RunMode, LogFile
from BALSAMIC.constants.cluster import (
    QOS,
    ClusterMailType,
)
from BALSAMIC.models.config import ConfigModel
from BALSAMIC.models.snakemake import SnakemakeExecutable
from BALSAMIC.utils.analysis import get_singularity_bind_paths
from BALSAMIC.utils.cli import createDir, get_snakefile
from BALSAMIC.utils.io import write_json, read_yaml, write_yaml
from BALSAMIC.utils.logging import add_file_logging

LOG = logging.getLogger(__name__)


@click.command("analysis", short_help="Run the analysis on a sample config-file")
@OPTION_CLUSTER_ACCOUNT
@OPTION_CLUSTER_MAIL
@OPTION_CLUSTER_MAIL_TYPE
@OPTION_CLUSTER_PROFILE
@OPTION_WORKFLOW_PROFILE
@OPTION_CLUSTER_QOS
@OPTION_CLUSTER_ENV
@OPTION_DRAGEN
@OPTION_FORCE_ALL
@OPTION_QUIET
@OPTION_RUN_ANALYSIS
@OPTION_RUN_MODE
@OPTION_RUN_INTERACTIVELY
@OPTION_SAMPLE_CONFIG
@OPTION_SNAKEFILE
@OPTION_SNAKEMAKE_OPT
@click.pass_context
def analysis(
    context: click.Context,
    snakefile: Path,
    sample_config: Path,
    run_mode: RunMode,
    cluster_env: Path,
    dragen: bool,
    cluster_profile: Path,
    workflow_profile: Path,
    run_analysis: bool,
    run_interactively: bool,
    qos: QOS,
    force_all: bool,
    snakemake_opt: List[str],
    account: str,
    mail_user: str,
    mail_type: ClusterMailType,
    quiet: bool,
):
    """Run BALSAMIC workflow on the provided sample's config file."""

    LOG.info(f"Initializing balsamic config model from config JSON: {sample_config}.")
    sample_config_path: Path = Path(sample_config).absolute()
    with open(sample_config_path, "r") as sample_fh:
        sample_config = json.load(sample_fh)

    config_model = ConfigModel.model_validate(sample_config)
    case_id = config_model.analysis.case_id

    log_file = Path(config_model.analysis.analysis_dir, case_id, LogFile.LOGNAME).as_posix()
    LOG.info(f"Setting BALSAMIC logfile path to: {log_file}.")
    add_file_logging(log_file, logger_name=__name__)

    LOG.info(f"Running BALSAMIC version {balsamic_version} -- RUN ANALYSIS")
    LOG.info(f"BALSAMIC started with log level {context.obj['log_level']}.")
    LOG.info(f"Using case config file: {sample_config_path}")
    LOG.info(f"Starting analysis on: {case_id}.")

    if run_mode == RunMode.CLUSTER and not run_analysis:
        LOG.info("Changing run-mode to local on dry-run")
        run_mode: RunMode = RunMode.LOCAL

    if run_mode == RunMode.CLUSTER and not account:
        LOG.info("An account is required for cluster run mode")
        raise click.Abort()

    # Create directories for results, logs, scripts and benchmark files
    result_path: Path = Path(config_model.analysis.result)
    log_path: Path = Path(config_model.analysis.log)
    script_path: Path = Path(config_model.analysis.script)
    benchmark_path: Path = Path(config_model.analysis.benchmark)

    LOG.info(f"Creating analysis and log directories.")
    analysis_directories_list = [result_path, log_path, script_path, benchmark_path]

    for analysis_sub_dir in analysis_directories_list:
        analysis_sub_dir.mkdir(exist_ok=True)

    if run_analysis:
        # if not dry run, and current existing log-dir is not empty, then create (new) log/script directory
        existing_log_files = os.listdir(log_path.as_posix())
        if existing_log_files:
            log_path = Path(createDir(log_path.as_posix(), []))
            script_path = Path(createDir(script_path.as_posix(), []))

    LOG.info(f"Updating config model with account: {account} and QOS: {qos}")
    config_model.qos = qos
    config_model.account = account
    config_model.analysis.log = log_path.as_posix()
    config_model.analysis.script = script_path.as_posix()

    config_model_dict: dict = config_model.model_dump(by_alias=True, exclude_none=True)
    LOG.info(f"Dumping updated config model to JSON: {sample_config_path}")
    write_json(json_obj=config_model_dict, path=sample_config_path)

    analysis_type = config_model.analysis.analysis_type
    analysis_workflow = config_model.analysis.analysis_workflow
    analysis_dir: Path = Path(config_model.analysis.analysis_dir)
    snakefile: Path = (
        snakefile if snakefile else get_snakefile(analysis_type, analysis_workflow)
    )
    LOG.info("Creating cluster profile yaml")
    cluster_yaml = read_yaml(Path(cluster_profile, "config.yaml").as_posix())

    # Update the placeholder fields
    cluster_yaml["default-resources"]["slurm_account"] = account
    cluster_yaml["default-resources"]["slurm_extra"] = f"--qos {qos} --error {log_path.as_posix()}/BALSAMIC.{case_id}.%x.%j.err --output {log_path.as_posix()}/BALSAMIC.{case_id}.%x.%j.stdout"

    write_yaml(cluster_yaml, Path(config_model.analysis.analysis_dir, case_id, "config.yaml").as_posix())

    LOG.info("Organizing snakemake run information")
    snakemake_executable: SnakemakeExecutable = SnakemakeExecutable(
        account=account,
        case_id=case_id,
        config_path=sample_config_path,
        dragen=dragen,
        force=force_all,
        log_dir=log_path.as_posix(),
        mail_type=mail_type,
        mail_user=mail_user,
        cluster_profile=Path(config_model.analysis.analysis_dir, case_id).as_posix(),
        workflow_profile=workflow_profile,
        qos=qos,
        quiet=quiet,
        run_analysis=run_analysis,
        run_mode=run_mode,
        script_dir=script_path.as_posix(),
        singularity_bind_paths=get_singularity_bind_paths(sample_config),
        snakefile=snakefile,
        snakemake_options=snakemake_opt,
        working_dir=Path(analysis_dir, case_id, "BALSAMIC_run"),
    )


    if not run_interactively:
        LOG.info(f"Creating sbatch-script to submit jobs.")
        # Create sbatch script here with snakemake run
        conda_env = os.environ.get("CONDA_DEFAULT_ENV", "")
        LOG.info(f"Setting conda environment in job-submission script to: {conda_env}")

        with open(f"{script_path.as_posix()}/BALSAMIC_snakemake_submit.sh", "w") as submit_file:
            sbatch_lines = ["#!/bin/bash -l",
                            f"#SBATCH --account={account}",
                            f"#SBATCH --job-name=BALSAMIC_snakemake_submit.{case_id}.%j",
                            f"#SBATCH --output={log_path}/BALSAMIC_snakemake_submit.{case_id}.%j.out",
                            f"#SBATCH --error={log_path}/BALSAMIC_snakemake_submit.{case_id}.%j.err",
                            "#SBATCH --ntasks=1",
                            "#SBATCH --mem=5G",
                            "#SBATCH --time=60:00:00",
                            f"#SBATCH --qos={qos}",
                            "#SBATCH --cpus-per-task=1"]

            if cluster_env:
                LOG.info(f"Setting cluster environment in job-submission script to: {cluster_env}")
                sbatch_lines.append(f"source {cluster_env}")

            sbatch_lines.extend([f"conda activate {conda_env}", f"{snakemake_executable.get_command()}"])
            submit_file.write("\n".join(sbatch_lines) + "\n")

        # Submit sbatch script to cluster
        subprocess.run(f"sbatch {script_path.as_posix()}/BALSAMIC_snakemake_submit.sh", shell=True)

    else:
        LOG.info(f"Starting {analysis_workflow} workflow...")
        subprocess.run(
            f"{sys.executable} -m {snakemake_executable.get_command()}",
            shell=True,
        )


