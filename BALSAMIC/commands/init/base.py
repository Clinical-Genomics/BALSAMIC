"""Balsamic init command."""
import json
import logging
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Union, List, Optional

import click

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
from BALSAMIC.constants.analysis import BIOINFO_TOOL_ENV, RunMode
from BALSAMIC.constants.cache import GenomeVersion, ContainerVersion, REFERENCE_FILES
from BALSAMIC.constants.cluster import ClusterMailType, QOS, ClusterProfile
from BALSAMIC.constants.paths import BALSAMIC_DIR
from BALSAMIC.models.cache import CacheConfig, ReferencesHg, ReferencesCanFam
from BALSAMIC.models.snakemake import Snakemake, SingularityBindPath
from BALSAMIC.utils.cli import get_snakefile
from BALSAMIC.utils.io import write_json, generate_workflow_graph

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
    context: click.Context,
    out_dir: str,
    container_version: ContainerVersion,
    cosmic_key: str,
    genome_version: GenomeVersion,
    snakefile: Path,
    run_mode: RunMode,
    cluster_config: Path,
    profile: ClusterProfile,
    qos: QOS,
    account: Optional[str],
    mail_user: Optional[str],
    mail_type: Optional[ClusterMailType],
    snakemake_opt: List[str],
    force_all: bool,
    run_analysis: bool,
    quiet: bool,
) -> None:
    """Validate inputs and download reference caches and containers."""
    LOG.info(f"BALSAMIC started with log level {context.obj['loglevel']}")

    if run_mode == RunMode.CLUSTER and not run_analysis:
        LOG.info("Changing run-mode to local on dry-run")
        run_mode: RunMode = RunMode.LOCAL

    if run_mode == RunMode.CLUSTER and not account:
        LOG.error("A cluster account is required for cluster run mode")
        raise click.Abort()

    if genome_version in [GenomeVersion.HG19, GenomeVersion.HG38] and not cosmic_key:
        LOG.error(
            f"No COSMIC authentication key specified. It is required when using {genome_version} reference"
        )
        raise click.Abort()

    out_dir: Path = Path(out_dir, balsamic_version).absolute()
    references_dir: Path = Path(out_dir, genome_version)
    genome_dir = Path(references_dir, "genome")
    variants_dir = Path(references_dir, "variants")
    vep_dir = Path(references_dir, "vep")
    containers_dir: Path = Path(out_dir, "containers")
    config_path: Path = Path(references_dir, "config.json")
    log_dir: Path = Path(references_dir, "logs")
    script_dir: Path = Path(references_dir, "scripts")
    for dir_path in [references_dir, log_dir, script_dir]:
        dir_path.mkdir(parents=True, exist_ok=True)

    references: Union[ReferencesHg, ReferencesCanFam] = REFERENCE_FILES[genome_version]
    cache_config: CacheConfig = CacheConfig(
        analysis={"case_id": f"reference.{genome_version}.v{balsamic_version}"},
        references_dir=references_dir.as_posix(),
        genome_dir=genome_dir.as_posix(),
        variants_dir=variants_dir.as_posix(),
        vep_dir=vep_dir.as_posix(),
        containers_dir=containers_dir.as_posix(),
        genome_version=genome_version,
        cosmic_key=cosmic_key,
        bioinfo_tools=BIOINFO_TOOL_ENV,
        containers=get_containers(container_version),
        references=references,
        references_date=datetime.now().strftime("%Y-%m-%d %H:%M"),
    )
    write_json(json.loads(cache_config.json(exclude_none=True)), config_path.as_posix())
    LOG.info(f"Reference workflow configured successfully ({config_path.as_posix()})")

    snakefile: Path = (
        snakefile if snakefile else get_snakefile("generate_ref", "balsamic")
    )

    generate_workflow_graph(
        config_path=config_path,
        directory_path=references_dir,
        snakefile=snakefile,
        title="reference",
    )

    LOG.info("Starting reference generation workflow...")
    snakemake: Snakemake = Snakemake(
        account=account,
        case_id=cache_config.analysis.case_id,
        cluster_config_path=cluster_config,
        config_path=config_path,
        force=force_all,
        log_dir=log_dir,
        mail_type=mail_type,
        mail_user=mail_user,
        profile=profile,
        qos=qos,
        quiet=quiet,
        result_dir=references_dir,
        run_analysis=run_analysis,
        run_mode=run_mode,
        script_dir=script_dir,
        snakemake_options=snakemake_opt,
        singularity=True,
        singularity_bind_paths=[
            SingularityBindPath(source=BALSAMIC_DIR, destination=BALSAMIC_DIR),
            SingularityBindPath(source=references_dir, destination=references_dir),
        ],
        snakefile=snakefile,
        working_dir=references_dir,
    )
    subprocess.run(
        f"{sys.executable} -m {snakemake.get_snakemake_command()}", shell=True
    )
