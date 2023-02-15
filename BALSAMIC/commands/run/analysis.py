import sys
import os
import logging
import subprocess
import json
import click

from pathlib import Path

from BALSAMIC.utils.cli import (
    createDir,
    get_schedulerpy,
    get_snakefile,
    SnakeMake,
    get_config,
    job_id_dump_to_yaml,
    get_input_files_path,
)
from BALSAMIC.constants.common import BALSAMIC_SCRIPTS
from BALSAMIC.constants.workflow_params import VCF_DICT

LOG = logging.getLogger(__name__)


@click.command(
    "analysis", short_help="Run the analysis on a provided sample config-file"
)
@click.option(
    "-S",
    "--snake-file",
    type=click.Path(),
    show_default=True,
    help=(
        "Input for a custom snakefile. WARNING: "
        "This is for internal testing, and should "
        "not be used. Providing a snakefile supersedes"
        "analysis_type option."
    ),
)
@click.option(
    "-s",
    "--sample-config",
    required=True,
    type=click.Path(),
    help="Sample json config file.",
)
@click.option(
    "--run-mode",
    show_default=True,
    default="cluster",
    type=click.Choice(["local", "cluster"]),
    help=(
        "Run mode to use. By default SLURM will be used to "
        "run the analysis. But local runner also available "
        "for local computing"
    ),
)
@click.option(
    "-c",
    "--cluster-config",
    show_default=True,
    default=get_config("cluster"),
    type=click.Path(),
    help="cluster config json file. (eg- SLURM, QSUB)",
)
@click.option(
    "--dragen", is_flag=True, default=False, help="Enable dragen variant caller"
)
@click.option(
    "-p",
    "--profile",
    default="slurm",
    type=click.Choice(["slurm", "qsub"]),
    help="cluster profile to submit jobs",
)
@click.option(
    "--benchmark",
    default=False,
    is_flag=True,
    help="Profile slurm jobs using the value of this option. Make sure you have slurm profiler enabled in your HPC.",
)
@click.option(
    "-r",
    "--run-analysis",
    show_default=True,
    default=False,
    is_flag=True,
    help=(
        "By default balsamic run_analysis will run in "
        "dry run mode. Raise thise flag to make the "
        "actual analysis"
    ),
)
@click.option(
    "--qos",
    type=click.Choice(["low", "normal", "high", "express"]),
    show_default=True,
    default="low",
    help="QOS for sbatch jobs. Passed to " + get_schedulerpy(),
)
@click.option(
    "-f",
    "--force-all",
    show_default=True,
    default=False,
    is_flag=True,
    help="Force run all analysis. This is same as snakemake --forceall",
)
@click.option(
    "--snakemake-opt", multiple=True, help="Pass these options directly to snakemake"
)
@click.option(
    "--account",
    "--slurm-account",
    "--qsub-account",
    help="cluster account to run jobs, ie: slurm_account",
)
@click.option(
    "--mail-user", help="cluster mail user to send out email. e.g.: slurm_mail_user"
)
@click.option(
    "-q",
    "--quiet",
    default=False,
    is_flag=True,
    help="Instruct snakemake to be quiet! No output will be printed",
)
@click.option(
    "--mail-type",
    type=click.Choice(
        [
            "NONE",
            "BEGIN",
            "END",
            "FAIL",
            "REQUEUE",
            "ALL",
            "TIME_LIMIT",
        ]
    ),
    help=(
        "cluster mail type to send out email. This will "
        "be applied to all jobs and override snakemake settings."
    ),
)
@click.option(
    "--disable-variant-caller",
    help=(
        f"Run workflow with selected variant caller(s) disable."
        f"Use comma to remove multiple variant callers. Valid "
        f"values are: {list(VCF_DICT.keys())}"
    ),
)
@click.pass_context
def analysis(
    context,
    snake_file,
    sample_config,
    run_mode,
    cluster_config,
    run_analysis,
    force_all,
    snakemake_opt,
    mail_type,
    mail_user,
    account,
    qos,
    profile,
    disable_variant_caller,
    quiet,
    dragen,
    benchmark,
):
    """
    Runs BALSAMIC workflow on the provided sample's config file
    """
    LOG.info(f"BALSAMIC started with log level {context.obj['loglevel']}.")

    if run_mode == "cluster" and not run_analysis:
        LOG.info("Changing run-mode to local on dry-run")
        run_mode = "local"

    if run_mode == "cluster" and not account:
        LOG.info(
            "slurm-account, qsub-account, or account is required for slurm run mode"
        )
        raise click.Abort()

    sample_config_path = os.path.abspath(sample_config)

    with open(sample_config, "r") as sample_fh:
        sample_config = json.load(sample_fh)

    logpath = sample_config["analysis"]["log"]
    scriptpath = sample_config["analysis"]["script"]
    resultpath = sample_config["analysis"]["result"]
    benchmarkpath = sample_config["analysis"]["benchmark"]
    case_name = sample_config["analysis"]["case_id"]

    if run_analysis:
        # if not dry run, then create (new) log/script directory
        for dirpath, dirnames, files in os.walk(logpath):
            if files:
                logpath = createDir(logpath, [])
                scriptpath = createDir(scriptpath, [])
                sample_config["analysis"]["benchmark"] = createDir(benchmarkpath, [])

    # Create result directory
    os.makedirs(resultpath, exist_ok=True)

    if not os.path.exists(logpath):
        os.makedirs(logpath, exist_ok=True)
        os.makedirs(scriptpath, exist_ok=True)
        os.makedirs(benchmarkpath, exist_ok=True)

    analysis_type = sample_config["analysis"]["analysis_type"]
    analysis_workflow = sample_config["analysis"]["analysis_workflow"]
    reference_genome = sample_config["reference"]["reference_genome"]

    # Singularity bind path
    bind_path = list()
    bind_path.extend(get_input_files_path(sample_config["analysis"]["fastq_path"]))
    bind_path.append(str(Path(__file__).parents[2] / "assets"))
    bind_path.append(os.path.commonpath(sample_config["reference"].values()))
    if "panel" in sample_config:
        bind_path.append(sample_config.get("panel").get("capture_kit"))
    if "background_variants" in sample_config:
        bind_path.append(sample_config.get("background_variants"))
    if "pon_cnn" in sample_config:
        bind_path.append(sample_config.get("panel").get("pon_cnn"))
    bind_path.append(BALSAMIC_SCRIPTS)
    bind_path.append(sample_config["analysis"]["analysis_dir"])

    # Construct snakemake command to run workflow
    balsamic_run = SnakeMake()
    balsamic_run.case_name = case_name
    balsamic_run.working_dir = (
        Path(
            sample_config["analysis"]["analysis_dir"], case_name, "BALSAMIC_run"
        ).as_posix()
        + "/"
    )
    balsamic_run.snakefile = (
        snake_file
        if snake_file
        else get_snakefile(analysis_type, analysis_workflow, reference_genome)
    )
    balsamic_run.configfile = sample_config_path
    balsamic_run.run_mode = run_mode
    balsamic_run.cluster_config = cluster_config
    balsamic_run.scheduler = get_schedulerpy()
    balsamic_run.profile = profile
    balsamic_run.log_path = logpath
    balsamic_run.script_path = scriptpath
    balsamic_run.result_path = resultpath
    balsamic_run.qos = qos
    balsamic_run.account = account
    if mail_type:
        balsamic_run.mail_type = mail_type
    balsamic_run.mail_user = mail_user
    balsamic_run.forceall = force_all
    balsamic_run.run_analysis = run_analysis
    balsamic_run.quiet = quiet
    # Always use singularity
    balsamic_run.use_singularity = True
    balsamic_run.singularity_bind = bind_path
    balsamic_run.sm_opt = snakemake_opt
    if benchmark and profile == "slurm":
        balsamic_run.slurm_profiler = "task"

    if disable_variant_caller:
        balsamic_run.disable_variant_caller = disable_variant_caller

    if dragen:
        balsamic_run.dragen = dragen

    cmd = sys.executable + " -m  " + balsamic_run.build_cmd()
    subprocess.run(cmd, shell=True)

    if run_analysis and run_mode == "cluster":
        jobid_dump = os.path.join(
            logpath, sample_config["analysis"]["case_id"] + ".sacct"
        )
        jobid_yaml = os.path.join(resultpath, profile + "_jobids.yaml")
        job_id_dump_to_yaml(jobid_dump, jobid_yaml, case_name)
