"""Balsamic initialize command."""
import json
import logging
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import click
import graphviz
import snakemake
from graphviz import Source

from BALSAMIC import __version__ as balsamic_version

from BALSAMIC.commands.init.options import (
    OPTION_OUT_DIR,
    OPTION_CONTAINER_VERSION,
    OPTION_COSMIC_KEY,
    OPTION_SNAKEFILE,
    OPTION_CLUSTER_CONFIG,
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
    OPTION_CLUSTER_MAIL,
    OPTION_CLUSTER_MAIL_TYPE,
    OPTION_QUIET,
)
from BALSAMIC.constants.analysis import RunMode
from BALSAMIC.constants.cache import GenomeVersion
from BALSAMIC.constants.common import BIOINFO_TOOL_ENV
from BALSAMIC.constants.paths import BALSAMIC_DIR
from BALSAMIC.constants.references import REFERENCE_FILES
from BALSAMIC.models.cache_models import CacheConfigModel
from BALSAMIC.utils.cli import SnakeMake, get_schedulerpy, get_snakefile, CaptureStdout
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

    out_dir: Path = Path(out_dir, balsamic_version)
    references_dir: Path = Path(out_dir, genome_version)
    genome_dir = Path(references_dir, "genome")
    variants_dir = Path(references_dir, "variants")
    vep_dir = Path(references_dir, "vep")
    containers_dir: Path = Path(out_dir, "containers")
    config_path: Path = Path(references_dir, "config.json")
    log_dir: Path = Path(references_dir, "logs")
    script_dir: Path = Path(references_dir, "scripts")
    for dir_path in [
        references_dir,
        genome_dir,
        variants_dir,
        vep_dir,
        containers_dir,
        log_dir,
        script_dir,
    ]:
        dir_path.mkdir(parents=True, exist_ok=True)

    cache_config: CacheConfigModel = CacheConfigModel(
        analysis={
            "case_id": "reference" + "." + genome_version + ".v" + balsamic_version
        },
        references_dir=references_dir.as_posix(),
        genome_dir=genome_dir.as_posix(),
        variants_dir=variants_dir.as_posix(),
        vep_dir=vep_dir.as_posix(),
        containers_dir=containers_dir.as_posix(),
        genome_version=genome_version,
        cosmic_key=cosmic_key,
        bioinfo_tools=BIOINFO_TOOL_ENV,
        containers=get_containers(container_version),
        references=REFERENCE_FILES[genome_version],
        references_date=datetime.now().strftime("%Y-%m-%d %H:%M"),
    )
    write_json(json.loads(cache_config.json(exclude_none=True)), config_path.as_posix())
    LOG.info(f"Reference workflow configured successfully ({config_path.as_posix()})")

    snakefile: Path = (
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
    graph_title: str = "_".join(["BALSAMIC", balsamic_version, "reference"])
    graph_dot: str = "".join(graph_dot).replace(
        "snakemake_dag {", 'BALSAMIC { label="' + graph_title + '";labelloc="t";'
    )
    graph: Source = graphviz.Source(
        graph_dot,
        directory=Path(references_dir).as_posix(),
        filename="reference_graph",
        format="pdf",
        engine="dot",
    )
    try:
        graph_pdf: Path = Path(graph.render())
        LOG.info(
            f"Reference workflow graph generated successfully ({graph_pdf.as_posix()}) "
        )
    except Exception:
        LOG.error("Reference workflow graph generation failed")
        raise click.Abort()

    LOG.info("Starting reference generation workflow...")
    balsamic_run: SnakeMake = SnakeMake()
    balsamic_run.account = account
    balsamic_run.case_name = cache_config.analysis.case_id
    balsamic_run.cluster_config = cluster_config
    balsamic_run.configfile = config_path.as_posix()
    balsamic_run.forceall = force_all
    balsamic_run.log_path = log_dir.as_posix()
    balsamic_run.mail_type = mail_type
    balsamic_run.mail_user = mail_user
    balsamic_run.profile = profile
    balsamic_run.qos = qos
    balsamic_run.quiet = quiet
    balsamic_run.result_path = references_dir.as_posix()
    balsamic_run.run_analysis = run_analysis
    balsamic_run.run_mode = run_mode
    balsamic_run.scheduler = get_schedulerpy()
    balsamic_run.script_path = script_dir.as_posix()
    balsamic_run.sm_opt = snakemake_opt
    balsamic_run.snakefile = snakefile
    balsamic_run.use_singularity = True
    balsamic_run.singularity_bind = [
        references_dir.as_posix(),
        BALSAMIC_DIR,
    ]
    balsamic_run.working_dir = references_dir

    cmd: str = sys.executable + " -m " + balsamic_run.build_cmd()
    subprocess.run(cmd, shell=True)
