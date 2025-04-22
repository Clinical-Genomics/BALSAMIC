"""Balsamic run analysis CLI."""
import json
import logging
import os
import re
import subprocess
import textwrap
import sys
from pathlib import Path
from typing import List

import click

from BALSAMIC import __version__ as balsamic_version

from BALSAMIC.commands.options import (
    OPTION_CLUSTER_ACCOUNT,
    OPTION_WORKFLOW_PROFILE,
    OPTION_CLUSTER_PROFILE,
    OPTION_MAX_RUN_HOURS,
    OPTION_CLUSTER_QOS,
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
)
from BALSAMIC.models.config import ConfigModel
from BALSAMIC.models.sbatchsubmitter import SbatchSubmitter
from BALSAMIC.models.snakemake import SnakemakeExecutable
from BALSAMIC.utils.analysis import get_singularity_bind_paths
from BALSAMIC.utils.cli import createDir, get_snakefile
from BALSAMIC.utils.io import write_json
from BALSAMIC.utils.logging import add_file_logging

LOG = logging.getLogger(__name__)


@click.command("analysis", short_help="Run the analysis on a sample config-file")
@OPTION_CLUSTER_ACCOUNT
@OPTION_CLUSTER_PROFILE
@OPTION_MAX_RUN_HOURS
@OPTION_WORKFLOW_PROFILE
@OPTION_CLUSTER_QOS
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
    dragen: bool,
    cluster_profile: Path,
    max_run_hours: int,
    workflow_profile: Path,
    run_analysis: bool,
    run_interactively: bool,
    qos: QOS,
    force_all: bool,
    snakemake_opt: List[str],
    account: str,
    quiet: bool,
):
    """Run BALSAMIC workflow on the provided sample's config file."""

    LOG.info(f"Initializing balsamic config model from config JSON: {sample_config}.")
    sample_config_path: Path = Path(sample_config).absolute()
    with open(sample_config_path, "r") as sample_fh:
        sample_config = json.load(sample_fh)

    config_model = ConfigModel.model_validate(sample_config)
    case_id = config_model.analysis.case_id

    log_file = Path(
        config_model.analysis.analysis_dir, case_id, LogFile.LOGNAME
    ).as_posix()
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

    LOG.info(f"Updating config model with account: {account}, QOS: {qos}")
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

    LOG.info("Organizing snakemake run information")
    snakemake_executable: SnakemakeExecutable = SnakemakeExecutable(
        account=account,
        case_id=case_id,
        config_path=sample_config_path,
        dragen=dragen,
        force=force_all,
        log_dir=log_path.as_posix(),
        cluster_profile=cluster_profile,
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
        LOG.info(f"Submitting {analysis_workflow} workflow to cluster.")
        submitter = SbatchSubmitter(
            case_id=case_id,
            script_path=Path(script_path),
            result_path=Path(result_path),
            log_path=Path(log_path),
            account=account,
            qos=qos,
            max_run_hours=max_run_hours,
            snakemake_executable=snakemake_executable,
            logger=LOG,
        )
        submitter.create_sbatch_script()
        job_id = submitter.submit_job()
        if job_id:
            submitter.write_job_id_yaml(job_id)

    else:
        LOG.info(f"Starting {analysis_workflow} workflow interactively.")
        subprocess.run(
            f"{sys.executable} -m {snakemake_executable.get_command()}",
            shell=True,
            check=True,
        )
