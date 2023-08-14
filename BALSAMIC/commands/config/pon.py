import os
import json
import logging
from pathlib import Path

import click

from BALSAMIC import __version__ as balsamic_version
from BALSAMIC.constants.paths import CONTAINERS_DIR
from BALSAMIC.utils.cli import (
    generate_graph,
    get_bioinfo_tools_version,
    get_pon_sample_list,
    get_analysis_fastq_files_directory,
)
from BALSAMIC.utils.io import write_json
from BALSAMIC.models.analysis import PonBalsamicConfigModel

from BALSAMIC.constants.analysis import BIOINFO_TOOL_ENV

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
    help="Path directing to list of PON fastq samples. NOTE: Fastq-files in directory requires this structure:"
    " X_X_X_[sampleID]_XXXXXX_R_2.fastq.gz",
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
    with open(reference_config, "r") as config_file:
        reference_dict = json.load(config_file)

    fastq_path: str = get_analysis_fastq_files_directory(
        case_dir=Path(analysis_dir, case_id).as_posix(), fastq_path=fastq_path
    )

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

    logpath = config_collection_dict["analysis"]["log"]
    scriptpath = config_collection_dict["analysis"]["script"]
    resultpath = config_collection_dict["analysis"]["result"]
    benchmarkpath = config_collection_dict["analysis"]["benchmark"]

    # Create result directory
    os.makedirs(resultpath, exist_ok=True)

    if not os.path.exists(logpath):
        os.makedirs(logpath, exist_ok=True)
        os.makedirs(scriptpath, exist_ok=True)
        os.makedirs(benchmarkpath, exist_ok=True)

    config_path = Path(analysis_dir, case_id, case_id + "_PON.json").as_posix()
    write_json(json_obj=config_collection_dict, path=config_path)
    LOG.info(f"PON config file saved successfully - {config_path}")

    try:
        generate_graph(config_collection_dict, config_path)
        LOG.info(f"BALSAMIC PON workflow has been configured successfully!")
    except ValueError:
        LOG.error(
            f'BALSAMIC PON dag graph generation failed - {config_collection_dict["analysis"]["dag"]}'
        )
        raise click.Abort()
