"""Balsamic panel of normals config case CLI."""
import json
import logging
import os
from pathlib import Path

import click
from BALSAMIC.commands.options import OPTION_GENOME_VERSION

from BALSAMIC import __version__ as balsamic_version
from BALSAMIC.commands.config.options import (
    OPTION_CASE_ID,
    OPTION_UMI,
    OPTION_UMI_TRIM_LENGTH,
    OPTION_QUALITY_TRIM,
    OPTION_ADAPTER_TRIM,
    OPTION_PANEL_BED,
    OPTION_BALSAMIC_CACHE,
    OPTION_ANALYSIS_DIR,
    OPTION_FASTQ_PATH,
    OPTION_PON_VERSION,
)
from BALSAMIC.constants.analysis import BIOINFO_TOOL_ENV
from BALSAMIC.constants.cache import GenomeVersion
from BALSAMIC.constants.paths import CONTAINERS_DIR
from BALSAMIC.models.analysis import ConfigModel
from BALSAMIC.utils.cli import (
    generate_graph,
    get_bioinfo_tools_version,
    get_pon_sample_list,
    get_analysis_fastq_files_directory,
)
from BALSAMIC.utils.io import write_json

LOG = logging.getLogger(__name__)


@click.command("pon", short_help="Create a sample config file for PON analysis")
@OPTION_ADAPTER_TRIM
@OPTION_ANALYSIS_DIR
@OPTION_BALSAMIC_CACHE
@OPTION_CASE_ID
@OPTION_FASTQ_PATH
@OPTION_GENOME_VERSION
@OPTION_PANEL_BED
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
    case_id: str,
    fastq_path: Path,
    genome_version: GenomeVersion,
    panel_bed: Path,
    quality_trim: bool,
    umi: bool,
    umi_trim_length: bool,
    version: str,
):
    reference_config = os.path.join(
        balsamic_cache, balsamic_version, genome_version, "reference.json"
    )
    with open(reference_config, "r") as config_file:
        reference_dict = json.load(config_file)

    fastq_path: str = get_analysis_fastq_files_directory(
        case_dir=Path(analysis_dir, case_id).as_posix(), fastq_path=fastq_path
    )

    config_collection_dict = ConfigModel(
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
            "analysis_type": "pon",
            "pon_version": version,
            "analysis_workflow": "balsamic",
            "sequencing_type": "targeted" if panel_bed else "wgs",
        },
        samples=get_pon_sample_list(fastq_path),
        reference=reference_dict,
        singularity={
            "image": os.path.join(balsamic_cache, balsamic_version, "containers")
        },
        bioinfo_tools=BIOINFO_TOOL_ENV,
        bioinfo_tools_version=get_bioinfo_tools_version(
            bioinfo_tools=BIOINFO_TOOL_ENV,
            container_conda_env_path=CONTAINERS_DIR,
        ),
        panel={"capture_kit": panel_bed} if panel_bed else None,
    ).dict(by_alias=True, exclude_none=True)
    LOG.info("PON config model instantiated successfully")

    result_path: Path = Path(config_collection_dict["analysis"]["result"])
    log_path: Path = Path(config_collection_dict["analysis"]["log"])
    script_path: Path = Path(config_collection_dict["analysis"]["script"])
    benchmark_path: Path = Path(config_collection_dict["analysis"]["benchmark"])

    # Create directories for results, logs, scripts and benchmark files
    analysis_directories_list = [result_path, log_path, script_path, benchmark_path]

    for analysis_sub_dir in analysis_directories_list:
        analysis_sub_dir.mkdir(exist_ok=True)

    config_path = Path(analysis_dir, case_id, case_id + "_PON.json").as_posix()
    write_json(json_obj=config_collection_dict, path=config_path)
    LOG.info(f"PON config file saved successfully - {config_path}")

    generate_graph(config_collection_dict, config_path)
    LOG.info(f"BALSAMIC PON workflow has been configured successfully!")
