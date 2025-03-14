"""Balsamic run analysis CLI."""
import json
import logging
import os
import subprocess
import sys
from pathlib import Path
from typing import List

import click

from BALSAMIC.commands.options import (
    OPTION_BENCHMARK,
    OPTION_CLUSTER_ACCOUNT,
    OPTION_CLUSTER_CONFIG,
    OPTION_CLUSTER_MAIL,
    OPTION_CLUSTER_MAIL_TYPE,
    OPTION_CLUSTER_PROFILE,
    OPTION_CLUSTER_QOS,
    OPTION_DISABLE_VARIANT_CALLER,
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
from BALSAMIC.constants.analysis import RunMode
from BALSAMIC.constants.cluster import (
    QOS,
    ClusterConfigType,
    ClusterMailType,
    ClusterProfile,
)
from BALSAMIC.models.config import ConfigModel
from BALSAMIC.models.snakemake import SnakemakeExecutable
from BALSAMIC.utils.analysis import get_singularity_bind_paths
from BALSAMIC.utils.cli import createDir, get_config_path, get_snakefile
from BALSAMIC.utils.io import write_sacct_to_yaml

LOG = logging.getLogger(__name__)


@click.command("analysis", short_help="Run the analysis on a sample config-file")
@OPTION_BENCHMARK
@OPTION_CLUSTER_ACCOUNT
@OPTION_CLUSTER_CONFIG
@OPTION_CLUSTER_MAIL
@OPTION_CLUSTER_MAIL_TYPE
@OPTION_CLUSTER_PROFILE
@OPTION_CLUSTER_QOS
@OPTION_DISABLE_VARIANT_CALLER
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
    cluster_config: Path,
    benchmark: bool,
    dragen: bool,
    profile: ClusterProfile,
    run_analysis: bool,
    run_interactively: bool,
    qos: QOS,
    force_all: bool,
    snakemake_opt: List[str],
    account: str,
    mail_user: str,
    mail_type: ClusterMailType,
    quiet: bool,
    disable_variant_caller: str,
):
    """Run BALSAMIC workflow on the provided sample's config file."""
    LOG.info(f"BALSAMIC started with log level {context.obj['log_level']}.")

    if run_mode == RunMode.CLUSTER and not run_analysis:
        LOG.info("Changing run-mode to local on dry-run")
        run_mode: RunMode = RunMode.LOCAL

    if run_mode == RunMode.CLUSTER and not account:
        LOG.info("An account is required for cluster run mode")
        raise click.Abort()

    sample_config_path: Path = Path(sample_config).absolute()
    with open(sample_config_path, "r") as sample_fh:
        sample_config = json.load(sample_fh)

    # Initialize balsamic model to run validation tests
    config_model = ConfigModel.model_validate(sample_config)

    case_name = config_model.analysis.case_id

    # Create directories for results, logs, scripts and benchmark files
    result_path: Path = Path(config_model.analysis.result)
    log_path: Path = Path(config_model.analysis.log)
    script_path: Path = Path(config_model.analysis.script)
    benchmark_path: Path = Path(config_model.analysis.benchmark)

    analysis_directories_list = [result_path, log_path, script_path, benchmark_path]

    for analysis_sub_dir in analysis_directories_list:
        analysis_sub_dir.mkdir(exist_ok=True)

    if run_analysis:
        # if not dry run, and current existing log-dir is not empty, then create (new) log/script directory
        existing_log_files = os.listdir(log_path.as_posix())
        if existing_log_files:
            log_path = Path(createDir(log_path.as_posix(), []))
            script_path = Path(createDir(script_path.as_posix(), []))

    for analysis_sub_dir in analysis_directories_list:
        analysis_sub_dir.mkdir(exist_ok=True)

    analysis_type = config_model.analysis.analysis_type
    analysis_workflow = config_model.analysis.analysis_workflow
    analysis_dir: Path = Path(config_model.analysis.analysis_dir)
    snakefile: Path = (
        snakefile if snakefile else get_snakefile(analysis_type, analysis_workflow)
    )

    LOG.info(f"Organizing snakemake run information")
    snakemake_executable: SnakemakeExecutable = SnakemakeExecutable(
        account=account,
        benchmark=benchmark,
        case_id=case_name,
        cluster_config_path=cluster_config
        if cluster_config
        else get_config_path(ClusterConfigType.ANALYSIS),
        config_path=sample_config_path,
        disable_variant_caller=disable_variant_caller,
        dragen=dragen,
        force=force_all,
        log_dir=log_path.as_posix(),
        mail_type=mail_type,
        mail_user=mail_user,
        profile=profile,
        qos=qos,
        quiet=quiet,
        run_analysis=run_analysis,
        run_mode=run_mode,
        script_dir=script_path.as_posix(),
        singularity_bind_paths=get_singularity_bind_paths(sample_config),
        snakefile=snakefile,
        snakemake_options=snakemake_opt,
        working_dir=Path(analysis_dir, case_name, "BALSAMIC_run"),
    )

    if not run_interactively:
        LOG.info(f"Creating sbatch-script to submit jobs.")
        # Create sbatch script here with snakemake run
        with open(f"{script_path.as_posix()}/BALSAMIC_snakemake_submit.sh", "w") as submit_file:
            submit_file.write("#!/bin/bash -l")
            submit_file.write("#SBATCH --account=development")
            submit_file.write("#SBATCH --ntasks=1")
            submit_file.write("#SBATCH --mem=5G")
            submit_file.write("#SBATCH --time=20:00:00")
            submit_file.write("#SBATCH --qos=high")
            submit_file.write("#SBATCH --cpus-per-task=1")
            submit_file.write("source /home/proj/production/servers/resources/hasta.scilifelab.se/usestage.sh")
            submit_file.write("conda activate D_balsamic_DEV_MJ")
            submit_file.write(f"{snakemake_executable.get_command()}")
        # Submit sbatch script to cluster

    else:
        LOG.info(f"Starting {analysis_workflow} workflow...")
        subprocess.run(
            f"{sys.executable} -m {snakemake_executable.get_command()}",
            shell=True,
        )


