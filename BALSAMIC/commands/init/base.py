"""Balsamic init command."""
import json
import logging
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Union

import click

from BALSAMIC.commands.options import (
    OPTION_CACHE_VERSION,
    OPTION_CLUSTER_ACCOUNT,
    OPTION_WORKFLOW_PROFILE,
    OPTION_CLUSTER_MAIL,
    OPTION_CLUSTER_MAIL_TYPE,
    OPTION_CLUSTER_PROFILE,
    OPTION_CLUSTER_QOS,
    OPTION_COSMIC_KEY,
    OPTION_FORCE_ALL,
    OPTION_GENOME_VERSION,
    OPTION_OUT_DIR,
    OPTION_QUIET,
    OPTION_RUN_ANALYSIS,
    OPTION_RUN_MODE,
    OPTION_SNAKEFILE,
    OPTION_SNAKEMAKE_OPT,
)
from BALSAMIC.constants.analysis import BIOINFO_TOOL_ENV, RunMode
from BALSAMIC.constants.cache import REFERENCE_FILES, GenomeVersion
from BALSAMIC.constants.cluster import (
    QOS,
    ClusterConfigType,
    ClusterMailType,
)
from BALSAMIC.models.cache import CacheConfig, ReferencesCanFam, ReferencesHg
from BALSAMIC.models.snakemake import SnakemakeExecutable
from BALSAMIC.utils.analysis import get_cache_singularity_bind_paths
from BALSAMIC.utils.cache import get_containers
from BALSAMIC.utils.cli import get_config_path, get_snakefile
from BALSAMIC.utils.io import generate_workflow_graph, write_json

LOG = logging.getLogger(__name__)


@click.command(
    "init", short_help="Download singularity containers and build the reference cache"
)
@OPTION_OUT_DIR
@OPTION_CACHE_VERSION
@OPTION_CLUSTER_ACCOUNT
@OPTION_CLUSTER_MAIL
@OPTION_CLUSTER_MAIL_TYPE
@OPTION_CLUSTER_PROFILE
@OPTION_WORKFLOW_PROFILE
@OPTION_CLUSTER_QOS
@OPTION_COSMIC_KEY
@OPTION_FORCE_ALL
@OPTION_GENOME_VERSION
@OPTION_QUIET
@OPTION_RUN_ANALYSIS
@OPTION_RUN_MODE
@OPTION_SNAKEFILE
@OPTION_SNAKEMAKE_OPT
@click.pass_context
def initialize(
    context: click.Context,
    account: Optional[str],
    cache_version: str,
    cosmic_key: str,
    force_all: bool,
    genome_version: GenomeVersion,
    mail_type: Optional[ClusterMailType],
    mail_user: Optional[str],
    out_dir: str,
    cluster_profile: Path,
    workflow_profile: Path,
    qos: QOS,
    quiet: bool,
    run_analysis: bool,
    run_mode: RunMode,
    snakefile: Path,
    snakemake_opt: List[str],
) -> None:
    """Validate inputs and download reference caches and containers."""
    LOG.info(f"BALSAMIC started with log level {context.obj['log_level']}")

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

    out_dir: Path = Path(out_dir, cache_version).absolute()
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
        analysis={"case_id": f"reference.{genome_version}.{cache_version}"},
        references_dir=references_dir.as_posix(),
        genome_dir=genome_dir.as_posix(),
        variants_dir=variants_dir.as_posix(),
        vep_dir=vep_dir.as_posix(),
        containers_dir=containers_dir.as_posix(),
        genome_version=genome_version,
        cosmic_key=cosmic_key,
        bioinfo_tools=BIOINFO_TOOL_ENV,
        containers=get_containers(cache_version),
        references=references,
        references_date=datetime.now().strftime("%Y-%m-%d %H:%M"),
    )
    write_json(
        json_obj=json.loads(cache_config.model_dump_json(exclude_none=True)),
        path=config_path.as_posix(),
    )
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
    snakemake_executable: SnakemakeExecutable = SnakemakeExecutable(
        account=account,
        case_id=cache_config.analysis.case_id,
        config_path=config_path,
        force=force_all,
        log_dir=log_dir,
        mail_type=mail_type,
        mail_user=mail_user,
        profile=cluster_profile,
        workflow_profile=workflow_profile,
        qos=qos,
        quiet=quiet,
        run_analysis=run_analysis,
        run_mode=run_mode,
        script_dir=script_dir,
        snakemake_options=snakemake_opt,
        singularity_bind_paths=get_cache_singularity_bind_paths(cache_config),
        snakefile=snakefile,
        working_dir=references_dir,
    )
    subprocess.run(
        f"{sys.executable} -m {snakemake_executable.get_command()}",
        shell=True,
    )
