import os
import json
import logging
from pathlib import Path

import click

from BALSAMIC import __version__ as balsamic_version
from BALSAMIC.utils.cli import (get_sample_dict, get_panel_chrom,
                                get_bioinfo_tools_version,
                                create_fastq_symlink, generate_graph)
from BALSAMIC.utils.constants import (CONTAINERS_CONDA_ENV_PATH, VCF_DICT,
                                      BIOINFO_TOOL_ENV)
from BALSAMIC.utils.models import BalsamicConfigModel

LOG = logging.getLogger(__name__)


@click.command("case",
               short_help="Create a sample config file from input sample data")
@click.option("--case-id",
              required=True,
              help="Sample id that is used for reporting, \
              naming the analysis jobs, and analysis path")
@click.option("--umi/--no-umi",
              default=True,
              show_default=True,
              is_flag=True,
              help=("UMI processing steps for samples with UMI tags."
                    "For WGS cases, UMI is always disabled."))
@click.option("--umi-trim-length",
              default=5,
              show_default=True,
              type=int,
              help="Trim N bases from reads in fastq")
@click.option("--quality-trim/--no-quality-trim",
              default=True,
              show_default=True,
              is_flag=True,
              help="Trim low quality reads in fastq")
@click.option("--adapter-trim/--no-adapter-trim",
              default=True,
              show_default=True,
              is_flag=True,
              help="Trim adapters from reads in fastq")
@click.option("-p",
              "--panel-bed",
              type=click.Path(exists=True, resolve_path=True),
              required=False,
              help="Panel bed file for variant calling.")
@click.option("-b",
              "--background-variants",
              type=click.Path(exists=True, resolve_path=True),
              required=False,
              help="Background set of valid variants for UMI")
@click.option("--balsamic-cache",
              type=click.Path(exists=True, resolve_path=True),
              required=True,
              help="Path to BALSAMIC cache")
@click.option("--analysis-dir",
              type=click.Path(exists=True, resolve_path=True),
              required=True,
              help="Root analysis path to store analysis logs and results. \
                                     The final path will be analysis-dir/sample-id"
              )
@click.option("-t",
              "--tumor",
              type=click.Path(exists=True, resolve_path=True),
              required=True,
              multiple=True,
              help="Fastq files for tumor sample.")
@click.option("-n",
              "--normal",
              type=click.Path(exists=True, resolve_path=True),
              required=False,
              multiple=True,
              help="Fastq files for normal sample.")
@click.option("--umiworkflow/--no-umiworkflow",
              default=True,
              show_default=True,
              is_flag=True,
              help="Enable running UMI workflow")
@click.option("--tumor-sample-name", help="Tumor sample name")
@click.option("--normal-sample-name", help="Normal sample name")
@click.option("-g",
              "--genome-version",
              default="hg19",
              type=click.Choice(["hg19", "hg38"]),
              help=("Genome version to prepare reference. Path to genome"
                    "will be <outdir>/genome_version"))
@click.pass_context
def case_config(context, case_id, umi, umi_trim_length, adapter_trim,
                quality_trim, panel_bed, background_variants, analysis_dir,
                tumor, normal, umiworkflow, tumor_sample_name,
                normal_sample_name, genome_version, balsamic_cache):

    try:
        samples = get_sample_dict(
            tumor=tumor,
            normal=normal,
            tumor_sample_name=tumor_sample_name,
            normal_sample_name=normal_sample_name,
        )
    except AttributeError:
        LOG.error(
            f"File name is invalid, use convention [SAMPLE_ID]_R_[1,2].fastq.gz"
        )
        raise click.Abort()

    reference_config = os.path.join(balsamic_cache, "reference",
                                    balsamic_version, genome_version,
                                    "reference.json")
    with open(reference_config, 'r') as f:
        reference_dict = json.load(f)["reference"]

    config_collection_dict = BalsamicConfigModel(
        QC={
            "quality_trim": quality_trim,
            "adapter_trim": adapter_trim,
            "umi_trim": umi if panel_bed else False,
            "umi_trim_length": umi_trim_length,
        },
        analysis={
            "case_id": case_id,
            "analysis_dir": analysis_dir,
            "analysis_type": "paired" if normal else "single",
            "sequencing_type": "targeted" if panel_bed else "wgs",
            "umiworkflow": umiworkflow
        },
        reference=reference_dict,
        singularity=os.path.join(balsamic_cache, "containers"),
        background_variants=background_variants,
        samples=samples,
        vcf=VCF_DICT,
        bioinfo_tools=BIOINFO_TOOL_ENV,
        bioinfo_tools_version=get_bioinfo_tools_version(
            BIOINFO_TOOL_ENV, CONTAINERS_CONDA_ENV_PATH),
        panel={
            "capture_kit": panel_bed,
            "chrom": get_panel_chrom(panel_bed),
        } if panel_bed else None,
    ).dict(by_alias=True, exclude_none=True)
    LOG.info("Config file generated successfully")

    Path.mkdir(Path(config_collection_dict["analysis"]["fastq_path"]),
               parents=True,
               exist_ok=True)
    LOG.info("Directories created successfully")

    create_fastq_symlink(
        casefiles=(tumor + normal),
        symlink_dir=Path(config_collection_dict["analysis"]["fastq_path"]),
    )
    LOG.info(f"Symlinks generated successfully")

    config_path = Path(analysis_dir) / case_id / (case_id + ".json")
    with open(config_path, "w+") as fh:
        fh.write(json.dumps(config_collection_dict, indent=4))
    LOG.info(f"Config file saved successfully - {config_path}")

    try:
        generate_graph(config_collection_dict, config_path)
        LOG.info(f"BALSAMIC Workflow has been configured successfully!")
    except ValueError as e:
        LOG.error(
            f'BALSAMIC dag graph generation failed - {config_collection_dict["analysis"]["dag"]}',
        )
        raise click.Abort()
