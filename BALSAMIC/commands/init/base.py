import os
import sys
import re
import logging
import subprocess
from pathlib import Path

import click
import graphviz
import snakemake

from BALSAMIC.constants.common import (
    BIOINFO_TOOL_ENV,
    BALSAMIC_DOCKER_PATH,
    VALID_CONTAINER_CONDA_NAME,
)
from BALSAMIC.utils.cli import (
    write_json,
    CaptureStdout,
    get_snakefile,
    SnakeMake,
    get_config,
    get_schedulerpy,
    # job_id_dump_to_yaml,
)
from BALSAMIC import __version__ as balsamic_version

LOG = logging.getLogger(__name__)


@click.command(
    "init", short_help="Download matching version for container and build reference"
)
@click.option(
    "-o",
    "--outdir",
    "--out-dir",
    required=True,
    help=(
        "Output directory for ref files."
        "This path will be used as base path for files"
    ),
)
@click.option(
    "-v",
    "--container-version",
    show_default=True,
    default=balsamic_version,
    type=click.Choice(["develop", "master", balsamic_version]),
    help="Container for BALSAMIC version to download",
)
@click.option(
    "-f",
    "--force",
    show_default=True,
    default=False,
    is_flag=True,
    help="Force re-downloading all containers",
)
@click.option("-c", "--cosmic-key", required=False, help="cosmic db authentication key")
@click.option(
    "-s",
    "--snakefile",
    default=None,
    type=click.Path(),
    show_default=True,
    help="snakefile for reference generation",
)
@click.option(
    "-d",
    "--dagfile",
    default="generate_ref_worflow_graph",
    show_default=True,
    help="DAG file for overview",
)
@click.option(
    "-g",
    "--genome-version",
    default="hg19",
    type=click.Choice(["hg19", "hg38", "canfam3"]),
    help=(
        "Genome version to prepare reference. Path to genome"
        "will be <outdir>/genome_version"
    ),
)
@click.option(
    "-r",
    "--run-analysis",
    show_default=True,
    default=False,
    is_flag=True,
    help=(
        "By default balsamic run_analysis will run in dry run mode."
        "Raise this flag to make the actual analysis"
    ),
)
@click.option(
    "--run-mode",
    show_default=True,
    default="cluster",
    type=click.Choice(["local", "cluster"]),
    help=(
        "Run mode to use. By default SLURM will be used to generate the balsamic_cache"
        "Alternatively, option for local computing"
    ),
)
@click.option(
    "--cluster-config",
    show_default=True,
    default=get_config("reference_cluster"),
    type=click.Path(),
    help="cluster config json file. (eg- SLURM, QSUB)",
)
@click.option(
    "-p",
    "--profile",
    default="slurm",
    type=click.Choice(["slurm", "qsub"]),
    help="cluster profile to submit jobs",
)
@click.option(
    "--account",
    "--slurm-account",
    "--qsub-account",
    help="cluster account to run jobs, ie: slurm_account",
)
@click.option(
    "--qos",
    type=click.Choice(["low", "normal", "high", "express"]),
    show_default=True,
    default="low",
    help="QOS for sbatch jobs. Passed to " + get_schedulerpy(),
)
@click.option(
    "--mail-user", help="cluster mail user to send out email. e.g.: slurm_mail_user"
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
    "-q",
    "--quiet",
    default=False,
    is_flag=True,
    help="Instruct snakemake to be quiet! No output will be printed",
)
@click.pass_context
def initialize(
    context,
    outdir,
    container_version,
    force,
    cosmic_key,
    snakefile,
    dagfile,
    genome_version,
    run_analysis,
    run_mode,
    cluster_config,
    account,
    qos,
    profile,
    mail_user,
    mail_type,
    force_all,
    quiet,
    snakemake_opt,
):
    """
    Initialize various resources after first installation.
    - Pull container(s) for BALSAMIC according to matching version
    - Download and build a reference
    """
    config_dict = dict()
    config_dict["singularity"] = dict()

    LOG.info("BALSAMIC started with log level %s" % context.obj["loglevel"])

    if run_mode == "cluster" and not run_analysis:
        LOG.info("Changing run-mode to local on dry-run")
        run_mode = "local"

    if run_mode == "cluster" and not account:
        LOG.info(
            "slurm-account, qsub-account, or account is required for slurm run mode"
        )
        raise click.Abort()

    if genome_version in ["hg38", "hg19"] and not cosmic_key:
        LOG.error("cosmic db authentication key required with hg38 and hg19")
        raise click.Abort()

    # resolve outdir to absolute path
    outdir = Path(outdir).resolve()
    container_outdir = Path(outdir, balsamic_version, "containers")
    Path(container_outdir).mkdir(parents=True, exist_ok=True)
    config_dict["singularity"]["image_path"] = container_outdir.as_posix()
    config_dict["singularity"]["containers"] = dict()

    pattern = re.compile(r"^(\d+\.)?(\d+\.)?(\*|\d+)$")
    if pattern.findall(container_version):
        docker_image_base_name = "release_v{}".format(container_version)
    else:
        docker_image_base_name = container_version

    for image_suffix in VALID_CONTAINER_CONDA_NAME:
        container_stub_url = "{}:{}-{}".format(
            BALSAMIC_DOCKER_PATH, docker_image_base_name, image_suffix
        )
        config_dict["singularity"]["containers"][image_suffix] = container_stub_url

    config_path = Path(__file__).parents[2] / "config"
    config_path = config_path.absolute()

    rule_directory = Path(__file__).parents[2]

    config_dict["bioinfo_tools"] = BIOINFO_TOOL_ENV
    config_dict["rule_directory"] = rule_directory.as_posix() + "/"

    reference_outdir = Path(outdir, balsamic_version, genome_version)
    Path(reference_outdir).mkdir(parents=True, exist_ok=True)
    config_json = Path(reference_outdir, "config.json").as_posix()
    dagfile_path = Path(reference_outdir, dagfile).as_posix()
    logpath = Path(reference_outdir, "logs")
    Path(logpath).mkdir(parents=True, exist_ok=True)
    scriptpath = Path(reference_outdir, "scripts")
    Path(scriptpath).mkdir(parents=True, exist_ok=True)

    config_dict["output"] = reference_outdir.as_posix()
    if cosmic_key:
        config_dict["cosmic_key"] = cosmic_key

    config_dict["genome_version"] = genome_version
    config_dict["analysis"] = {}
    config_dict["analysis"]["case_id"] = (
        "reference" + "." + genome_version + ".v" + balsamic_version
    )

    write_json(config_dict, config_json)
    LOG.info("Reference generation workflow configured successfully - %s" % config_json)

    snakefile = (
        snakefile
        if snakefile
        else get_snakefile("generate_ref", "balsamic", genome_version)
    )

    with CaptureStdout() as graph_dot:
        snakemake.snakemake(
            snakefile=snakefile,
            dryrun=True,
            configfiles=[config_json],
            printrulegraph=True,
        )

    graph_title = "_".join(["BALSAMIC", balsamic_version, "Generate reference"])
    graph_dot = "".join(graph_dot).replace(
        "snakemake_dag {", 'BALSAMIC { label="' + graph_title + '";labelloc="t";'
    )
    graph_obj = graphviz.Source(
        graph_dot, filename=dagfile_path, format="pdf", engine="dot"
    )

    try:
        graph_pdf = graph_obj.render()
        LOG.info("Reference workflow graph generated successfully - %s " % graph_pdf)
    except Exception:
        LOG.error("Reference workflow graph generation failed")
        raise click.Abort()

    LOG.info("Reference generation workflow started")

    # Singularity bind path
    bind_path = list()
    bind_path.append(config_dict["output"])
    bind_path.append(config_dict["rule_directory"])

    # Construct snakemake command to run workflow
    balsamic_run = SnakeMake()
    balsamic_run.working_dir = config_dict["output"]
    balsamic_run.snakefile = snakefile
    balsamic_run.configfile = config_json
    balsamic_run.run_mode = run_mode
    balsamic_run.forceall = force_all
    balsamic_run.run_analysis = run_analysis
    balsamic_run.cluster_config = cluster_config
    balsamic_run.scheduler = get_schedulerpy()
    balsamic_run.profile = profile
    balsamic_run.account = account
    balsamic_run.qos = qos
    balsamic_run.log_path = logpath
    balsamic_run.script_path = scriptpath
    balsamic_run.result_path = reference_outdir
    balsamic_run.case_name = config_dict["analysis"]["case_id"]
    balsamic_run.quiet = quiet
    if mail_type:
        balsamic_run.mail_type = mail_type
    balsamic_run.mail_user = mail_user
    balsamic_run.sm_opt = snakemake_opt

    # Always use singularity
    balsamic_run.use_singularity = True
    balsamic_run.singularity_bind = bind_path

    cmd = sys.executable + " -m " + balsamic_run.build_cmd()
    subprocess.run(cmd, shell=True)
