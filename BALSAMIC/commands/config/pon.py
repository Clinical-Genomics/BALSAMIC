import os
import json
import logging
from pathlib import Path

import click

from BALSAMIC import __version__ as balsamic_version
from BALSAMIC.utils.cli import generate_graph, get_bioinfo_tools_version
from BALSAMIC.utils.io import write_json, remove_symlinks
from BALSAMIC.utils.models import PonBalsamicConfigModel

from BALSAMIC.constants.common import (
    CONTAINERS_CONDA_ENV_PATH,
    BIOINFO_TOOL_ENV,
)

LOG = logging.getLogger(__name__)


@click.command("pon", short_help="Create a sample config file for PON analysis")
@click.option("--case-id", required=True, help="Sample id used for reporting analysis")
@click.option(
    "--umi/--no-umi",
    default=True,
    show_default=True,
    is_flag=True,
    help=(
        "UMI processing steps for samples with UMI tags."
        "For WGS cases,by default UMI is disabled."
    ),
)
@click.option(
    "--umi-trim-length",
    default=5,
    show_default=True,
    type=int,
    help="Trimming first N bases from reads in fastq file",
)
@click.option(
    "--quality-trim/--no-quality-trim",
    default=True,
    show_default=True,
    is_flag=True,
    help="Trimming low quality reads in fastq file",
)
@click.option(
    "--adapter-trim/--no-adapter-trim",
    default=True,
    show_default=True,
    is_flag=True,
    help="Preprocess fastq reads by trimming adapters",
)
@click.option(
    "-p",
    "--panel-bed",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    help="Panel bed file for calculating target regions.",
)
@click.option(
    "--balsamic-cache",
    type=click.Path(exists=True, resolve_path=True),
    required=True,
    help="Path to BALSAMIC cache",
)
@click.option(
    "--analysis-dir",
    type=click.Path(exists=True, resolve_path=True),
    required=True,
    help="Root analysis path directory.",
)
@click.option(
    "--fastq-path",
    type=click.Path(exists=True, resolve_path=True),
    required=True,
    help="Path directing to list of PON fastq samples.",
)
@click.option(
    "-g",
    "--genome-version",
    default="hg19",
    type=click.Choice(["hg19"]),
    help=(
        "Genome version to prepare reference. Path to genome"
        "will be <outdir>/genome_version"
    ),
)
@click.option(
    "-v",
    "--version",
    default="v1",
    type=str,
    help="Version of the PON file to be generated",
)
@click.pass_context
def pon_config(
    context,
    case_id,
    analysis_dir,
    fastq_path,
    panel_bed,
    quality_trim,
    umi,
    umi_trim_length,
    adapter_trim,
    genome_version,
    balsamic_cache,
    version,
):
    reference_config = os.path.join(
        balsamic_cache, balsamic_version, genome_version, "reference.json"
    )
    with open(reference_config, "r") as f:
        reference_dict = json.load(f)["reference"]

    config_collection_dict = PonBalsamicConfigModel(
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
        reference=reference_dict,
        singularity=os.path.join(balsamic_cache, balsamic_version, "containers"),
        bioinfo_tools=BIOINFO_TOOL_ENV,
        bioinfo_tools_version=get_bioinfo_tools_version(
            BIOINFO_TOOL_ENV, CONTAINERS_CONDA_ENV_PATH
        ),
        panel={"capture_kit": panel_bed} if panel_bed else None,
    ).dict(by_alias=True, exclude_none=True)
    LOG.info("PON config file generated successfully")

    config_path = Path(analysis_dir) / case_id / (case_id + "_PON" + ".json")
    write_json(json_obj=config_collection_dict, path=config_path)
    LOG.info(f"PON config file saved successfully - {config_path}")
    remove_symlinks(directory=fastq_path, pattern="*.fastq.gz")

    try:
        generate_graph(config_collection_dict, config_path)
        LOG.info(f"BALSAMIC PON workflow has been configured successfully!")
    except ValueError:
        LOG.error(
            f'BALSAMIC PON dag graph generation failed - {config_collection_dict["analysis"]["dag"]}'
        )
        raise click.Abort()
