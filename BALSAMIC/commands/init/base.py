"""Balsamic initialize command."""
import logging
import subprocess
import sys
from pathlib import Path

import click
import graphviz
import snakemake
from BALSAMIC import __version__ as balsamic_version

from BALSAMIC.commands.init.options import (
    OPTION_OUT_DIR,
    OPTION_CONTAINER_VERSION,
    OPTION_COSMIC_KEY,
    OPTION_SNAKEFILE,
)
from BALSAMIC.commands.init.utils import get_containers
from BALSAMIC.commands.options import (
    OPTION_RUN_MODE,
    OPTION_CLUSTER_PROFILE,
    OPTION_CLUSTER_QOS,
    OPTION_CLUSTER_ACCOUNT,
    OPTION_SNAKEMAKE_OPT,
    OPTION_GENOME_VERSION,
    OPTION_FORCE_ALL,
    OPTION_RUN_ANALYSIS,
    OPTION_CLUSTER_CONFIG,
    OPTION_CLUSTER_MAIL,
    OPTION_CLUSTER_MAIL_TYPE,
    OPTION_QUIET,
)
from BALSAMIC.constants.analysis import RunMode, GenomeVersion
from BALSAMIC.constants.common import BIOINFO_TOOL_ENV, ConfigType
from BALSAMIC.utils.cli import (
    SnakeMake,
    get_schedulerpy,
    get_snakefile,
    CaptureStdout,
    get_config_path,
)
from BALSAMIC.utils.io import write_json

LOG = logging.getLogger(__name__)


@click.command(
    "init", short_help="Download singularity containers and build the reference cache"
)
@OPTION_OUT_DIR
@OPTION_CONTAINER_VERSION
@OPTION_COSMIC_KEY
@OPTION_GENOME_VERSION
@OPTION_SNAKEFILE
@OPTION_RUN_MODE
@OPTION_CLUSTER_CONFIG
@OPTION_CLUSTER_PROFILE
@OPTION_CLUSTER_QOS
@OPTION_CLUSTER_ACCOUNT
@OPTION_CLUSTER_MAIL
@OPTION_CLUSTER_MAIL_TYPE
@OPTION_SNAKEMAKE_OPT
@OPTION_FORCE_ALL
@OPTION_RUN_ANALYSIS
@OPTION_QUIET
@click.pass_context
def initialize(
    context,
    out_dir,
    container_version,
    cosmic_key,
    genome_version,
    snakefile,
    run_mode,
    cluster_config,
    profile,
    qos,
    account,
    mail_user,
    mail_type,
    snakemake_opt,
    force_all,
    run_analysis,
    quiet,
):
    """Reference cache and containers download."""
    LOG.info(f"BALSAMIC started with log level {context.obj['loglevel']}")

    if run_mode == RunMode.CLUSTER and not run_analysis:
        LOG.info("Changing run-mode to local on dry-run")
        run_mode = RunMode.LOCAL

    if run_mode == RunMode.CLUSTER and not account:
        LOG.error("A cluster account is required for cluster run mode")
        raise click.Abort()

    if genome_version in [GenomeVersion.HG19, GenomeVersion.HG38] and not cosmic_key:
        LOG.error(
            f"A COSMIC authentication key is required for {GenomeVersion.HG19}/{GenomeVersion.HG38}"
        )
        raise click.Abort()

    out_dir: Path = Path(out_dir).resolve()
    containers_dir: Path = Path(out_dir, balsamic_version, "containers")
    reference_dir: Path = Path(out_dir, balsamic_version, genome_version)
    config_path: Path = Path(reference_dir, "config.json")
    log_dir: Path = Path(reference_dir, "logs")
    script_dir: Path = Path(reference_dir, "scripts")
    for dir_path in [containers_dir, reference_dir, log_dir, script_dir]:
        dir_path.mkdir(parents=True, exist_ok=True)

    config_dict: dict = {
        "output": reference_dir.as_posix(),
        "bioinfo_tools": BIOINFO_TOOL_ENV,
        "genome_version": genome_version,
        "rule_directory": Path(__file__).parents[2].as_posix() + "/",
        "singularity": {
            "image_path": containers_dir.as_posix(),
            "containers": get_containers(container_version),
        },
    }
    config_dict.update({"cosmic_key": cosmic_key}) if cosmic_key else None
    write_json(config_dict, config_path.as_posix())
    LOG.info(f"Reference workflow configured successfully ({config_path.as_posix()})")

    snakefile = (
        snakefile
        if snakefile
        else get_snakefile("generate_ref", "balsamic", genome_version)
    )

    with CaptureStdout() as graph_dot:
        snakemake.snakemake(
            snakefile=snakefile,
            dryrun=True,
            configfiles=[config_path.as_posix()],
            printrulegraph=True,
        )
    graph_title = "_".join(["BALSAMIC", balsamic_version, "reference"])
    graph_dot = "".join(graph_dot).replace(
        "snakemake_dag {", 'BALSAMIC { label="' + graph_title + '";labelloc="t";'
    )
    graph_obj = graphviz.Source(
        graph_dot,
        directory=Path(reference_dir).as_posix(),
        filename="reference_graph",
        format="pdf",
        engine="dot",
    )
    try:
        graph_pdf = graph_obj.render()
        LOG.info(f"Reference workflow graph generated successfully ({graph_pdf}) ")
    except Exception:
        LOG.error("Reference workflow graph generation failed")
        raise click.Abort()

    LOG.info("Starting reference generation workflow...")
    balsamic_run = SnakeMake()
    balsamic_run.configfile = config_path.as_posix()
    balsamic_run.case_name = "reference." + genome_version + ".v" + balsamic_version
    balsamic_run.working_dir = config_dict["output"]
    balsamic_run.snakefile = snakefile
    balsamic_run.run_mode = run_mode
    balsamic_run.scheduler = get_schedulerpy()
    balsamic_run.cluster_config = (
        cluster_config
        if cluster_config
        else get_config_path(ConfigType.CLUSTER_REFERENCE)
    )
    balsamic_run.profile = profile
    balsamic_run.qos = qos
    balsamic_run.account = account
    balsamic_run.mail_user = mail_user
    balsamic_run.mail_type = mail_type
    balsamic_run.sm_opt = snakemake_opt
    balsamic_run.forceall = force_all
    balsamic_run.run_analysis = run_analysis
    balsamic_run.quiet = quiet
    balsamic_run.log_path = log_dir.as_posix()
    balsamic_run.script_path = script_dir.as_posix()
    balsamic_run.result_path = reference_dir.as_posix()
    balsamic_run.use_singularity = True
    balsamic_run.singularity_bind = [
        config_dict["output"],
        config_dict["rule_directory"],
    ]

    cmd = sys.executable + " -m " + balsamic_run.build_cmd()
    subprocess.run(cmd, shell=True)
