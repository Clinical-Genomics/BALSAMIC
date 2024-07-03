"""Balsamic panel of normals config case CLI."""
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict

import click

from BALSAMIC import __version__ as balsamic_version
from BALSAMIC.commands.options import (
    OPTION_ADAPTER_TRIM,
    OPTION_ANALYSIS_DIR,
    OPTION_BALSAMIC_CACHE,
    OPTION_CACHE_VERSION,
    OPTION_CASE_ID,
    OPTION_FASTQ_PATH,
    OPTION_GENOME_INTERVAL,
    OPTION_GENOME_VERSION,
    OPTION_SENTIEON_INSTALL_DIR,
    OPTION_SENTIEON_LICENSE,
    OPTION_PANEL_BED,
    OPTION_PON_VERSION,
    OPTION_PON_WORKFLOW,
    OPTION_QUALITY_TRIM,
    OPTION_UMI,
    OPTION_UMI_TRIM_LENGTH,
)
from BALSAMIC.constants.analysis import BIOINFO_TOOL_ENV, PONWorkflow
from BALSAMIC.constants.cache import GenomeVersion
from BALSAMIC.constants.constants import FileType
from BALSAMIC.constants.paths import (CONTAINERS_DIR, SENTIEON_DNASCOPE_MODEL, SENTIEON_TNSCOPE_MODEL)
from BALSAMIC.models.config import ConfigModel
from BALSAMIC.utils.cli import (
    generate_graph,
    get_analysis_fastq_files_directory,
    get_bioinfo_tools_version,
    get_pon_sample_list,
)
from BALSAMIC.utils.io import read_json, write_json
from BALSAMIC.utils.utils import get_absolute_paths_dict

LOG = logging.getLogger(__name__)


@click.command("pon", short_help="Create a sample config file for PON analysis")
@OPTION_ADAPTER_TRIM
@OPTION_ANALYSIS_DIR
@OPTION_BALSAMIC_CACHE
@OPTION_CACHE_VERSION
@OPTION_CASE_ID
@OPTION_FASTQ_PATH
@OPTION_GENOME_VERSION
@OPTION_GENOME_INTERVAL
@OPTION_SENTIEON_INSTALL_DIR
@OPTION_SENTIEON_LICENSE
@OPTION_PANEL_BED
@OPTION_PON_WORKFLOW
@OPTION_PON_VERSION
@OPTION_QUALITY_TRIM
@OPTION_UMI
@OPTION_UMI_TRIM_LENGTH
@click.pass_context
def pon_config(
    context: click.Context,
    adapter_trim: bool,
    analysis_dir: Path,
    balsamic_cache: Path,
    cache_version: str,
    case_id: str,
    fastq_path: Path,
    genome_version: GenomeVersion,
    genome_interval: Path,
    sentieon_install_dir: Path,
    sentieon_license: str,
    panel_bed: Path,
    pon_workflow: PONWorkflow,
    quality_trim: bool,
    umi: bool,
    umi_trim_length: bool,
    version: str,
):
    references_path: Path = Path(balsamic_cache, cache_version, genome_version)
    references: Dict[str, Path] = get_absolute_paths_dict(
        base_path=references_path,
        data=read_json(Path(references_path, f"reference.{FileType.JSON}").as_posix()),
    )

    if pon_workflow in [PONWorkflow.GENS_MALE, PONWorkflow.GENS_FEMALE]:
        if not genome_interval:
            raise click.BadParameter(
                "Argument: genome_interval is required for GENS PON creation."
            )
        references["genome_interval"] = genome_interval

    if pon_workflow == PONWorkflow.CNVKIT and not panel_bed:
        raise click.BadParameter(
            "Argument: panel_bed is required for CNVkit PON creation."
        )

    fastq_path: str = get_analysis_fastq_files_directory(
        case_dir=Path(analysis_dir, case_id).as_posix(), fastq_path=fastq_path
    )
    result_dir: Path = Path(analysis_dir, case_id, "analysis")
    log_dir: Path = Path(analysis_dir, case_id, "logs")
    script_dir: Path = Path(analysis_dir, case_id, "scripts")
    benchmark_dir: Path = Path(analysis_dir, case_id, "benchmarks")
    dag_path: Path = Path(
        analysis_dir, case_id, f"{case_id}_BALSAMIC_{balsamic_version}_graph.pdf"
    )
    for directory in [result_dir, log_dir, script_dir, benchmark_dir]:
        directory.mkdir(exist_ok=True)

    config_collection_dict = ConfigModel(
        sentieon={
            "sentieon_install_dir": sentieon_install_dir,
            "sentieon_license": sentieon_license,
            "sentieon_exec": Path(sentieon_install_dir, "bin", "sentieon").as_posix(),
            "dnascope_model": SENTIEON_DNASCOPE_MODEL,
            "tnscope_model": SENTIEON_TNSCOPE_MODEL,
        },
        QC={
            "adapter_trim": adapter_trim,
            "quality_trim": quality_trim,
            "umi_trim": umi if panel_bed else False,
            "umi_trim_length": umi_trim_length,
        },
        analysis={
            "case_id": case_id,
            "analysis_dir": analysis_dir,
            "fastq_path": fastq_path,
            "log": log_dir.as_posix(),
            "script": script_dir.as_posix(),
            "result": result_dir.as_posix(),
            "benchmark": benchmark_dir.as_posix(),
            "dag": dag_path.as_posix(),
            "analysis_type": "pon",
            "pon_workflow": pon_workflow,
            "pon_version": version,
            "analysis_workflow": "balsamic",
            "sequencing_type": "targeted" if panel_bed else "wgs",
            "config_creation_date": datetime.now().strftime("%Y-%m-%d %H:%M"),
        },
        samples=get_pon_sample_list(fastq_path),
        reference=references,
        singularity={
            "image": Path(balsamic_cache, cache_version, "containers").as_posix()
        },
        bioinfo_tools=BIOINFO_TOOL_ENV,
        bioinfo_tools_version=get_bioinfo_tools_version(
            bioinfo_tools=BIOINFO_TOOL_ENV,
            container_conda_env_path=CONTAINERS_DIR,
        ),
        panel={"capture_kit": panel_bed} if panel_bed else None,
    ).model_dump(by_alias=True, exclude_none=True)
    LOG.info("PON config model instantiated successfully")

    config_path = Path(analysis_dir, case_id, case_id + "_PON.json").as_posix()
    write_json(json_obj=config_collection_dict, path=config_path)
    LOG.info(f"PON config file saved successfully - {config_path}")

    generate_graph(config_collection_dict, config_path)
    LOG.info(f"BALSAMIC PON workflow has been configured successfully!")
