import os
import json
import logging
from pathlib import Path

import click

from BALSAMIC import __version__ as balsamic_version
from BALSAMIC.utils.cli import (
    get_sample_dict,
    get_panel_chrom,
    get_bioinfo_tools_version,
    create_fastq_symlink,
    generate_graph,
)
from BALSAMIC.constants.common import (
    CONTAINERS_CONDA_ENV_PATH,
    BIOINFO_TOOL_ENV,
)
from BALSAMIC.constants.workflow_params import VCF_DICT
from BALSAMIC.utils.models import BalsamicConfigModel

LOG = logging.getLogger(__name__)

@click.command(
    "qc",
    short_help="Create a sample config file from input sample data for QC analysis",
)
@click.option(
    "--case-id",
    required=True,
    help="Sample id that is used for reporting, \
         naming the analysis jobs, and analysis path",
)
@click.option(
    "--quality-trim/--no-quality-trim",
    default=True,
    show_default=True,
    is_flag=True,
    help="Trim low quality reads in fastq",
)
@click.option(
    "--adapter-trim/--no-adapter-trim",
    default=True,
    show_default=True,
    is_flag=True,
    help="Trim adapters from reads in fastq",
)
@click.option(
    "-p",
    "--panel-bed",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    help="Panel bed file for variant calling.",
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
    help="Root analysis path to store analysis logs and results. \
          The final path will be analysis-dir/sample-id",
)
@click.option(
    "-t",
    "--tumor",
    type=click.Path(exists=True, resolve_path=True),
    required=True,
    multiple=True,
    help="Fastq files for tumor sample.",
)
@click.option(
    "-n",
    "--normal",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    multiple=True,
    help="Fastq files for normal sample.",
)
@click.option(
    "--umiworkflow/--no-umiworkflow",
    default=False,
    show_default=True,
    is_flag=True,
    help="Enable running UMI workflow",
)
@click.option("--tumor-sample-name", help="Tumor sample name")
@click.option("--normal-sample-name", help="Normal sample name")
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
@click.pass_context
def qc_config(
    context,
    case_id,
    adapter_trim,
    quality_trim,
    panel_bed,
    umiworkflow,
    analysis_dir,
    tumor,
    normal,
    tumor_sample_name,
    normal_sample_name,
    genome_version,
    balsamic_cache,
):

    try:
        samples = get_sample_dict(
            tumor=tumor,
            normal=normal,
            tumor_sample_name=tumor_sample_name,
            normal_sample_name=normal_sample_name,
        )
    except AttributeError:
        LOG.error(f"File name is invalid, use convention [SAMPLE_ID]_R_[1,2].fastq.gz")
        raise click.Abort()

    reference_config = os.path.join(
        balsamic_cache, balsamic_version, genome_version, "reference.json"
    )
    with open(reference_config, "r") as f:
        reference_dict = json.load(f)["reference"]

    config_collection_dict = BalsamicConfigModel(
        QC={
            "quality_trim": quality_trim,
            "adapter_trim": adapter_trim,
        },
        analysis={
            "case_id": case_id,
            "analysis_dir": analysis_dir,
            "analysis_type": "qc",
            "sequencing_type": "targeted" if panel_bed else "wgs",
        },
        reference=reference_dict,
        singularity=os.path.join(balsamic_cache, balsamic_version, "containers"),
        samples=samples,
        vcf=VCF_DICT,
        bioinfo_tools=BIOINFO_TOOL_ENV,
        bioinfo_tools_version=get_bioinfo_tools_version(
            BIOINFO_TOOL_ENV, CONTAINERS_CONDA_ENV_PATH
        ),
        panel={
            "capture_kit": panel_bed,
            "chrom": get_panel_chrom(panel_bed),
        }
        if panel_bed
        else None,
        umiworkflow=False,
    ).dict(by_alias=True, exclude_none=True)
    LOG.info("QC config file generated successfully")

    Path.mkdir(
        Path(config_collection_dict["analysis"]["fastq_path"]),
        parents=True,
        exist_ok=True,
    )
    LOG.info("Directories created successfully")

    create_fastq_symlink(
        casefiles=(tumor + normal),
        symlink_dir=Path(config_collection_dict["analysis"]["fastq_path"]),
    )
    LOG.info(f"Symlinks generated successfully")

    config_path = Path(analysis_dir) / case_id / (case_id + "_QC.json")
    with open(config_path, "w+") as fh:
        fh.write(json.dumps(config_collection_dict, indent=4))
    LOG.info(f"QC config file saved successfully - {config_path}")

    try:
        generate_graph(config_collection_dict, config_path)
        LOG.info(f"BALSAMIC QC Workflow has been configured successfully!")
    except ValueError:
        LOG.error(
            f'BALSAMIC QC dag graph generation failed - {config_collection_dict["analysis"]["dag"]}',
        )
        raise click.Abort()
