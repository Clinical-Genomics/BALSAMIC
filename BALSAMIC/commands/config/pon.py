import os
import json
import logging
from pathlib import Path

import click

from BALSAMIC import __version__ as balsamic_version
from BALSAMIC.utils.cli import (pon_sample_dict,
                                create_fastq_symlink, generate_graph)
from BALSAMIC.utils.models import PonBalsamicConfigModel

LOG = logging.getLogger(__name__)

@click.command("pon",
               short_help="Create a sample config file for PON analysis")
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
@click.option("-n",
              "--normal",
              type=click.Path(exists=True, resolve_path=True),
              required=True,
              multiple=True,
              help="Fastq files for normal sample.")
@click.option("--normal-sample-name", help="Normal sample name")
@click.option("-g",
              "--genome-version",
              default="hg19",
              type=click.Choice(["hg19", "hg38"]),
              help=("Genome version to prepare reference. Path to genome"
                    "will be <outdir>/genome_version"))
@click.pass_context
def pon_config(context, case_id, umi, umi_trim_length, adapter_trim,
                quality_trim, panel_bed, analysis_dir, normal, 
		normal_sample_name, genome_version, balsamic_cache):

    try:
        samples = pon_sample_dict(
            normal=normal,
            normal_sample_name=normal_sample_name
        )
    except AttributeError:
        LOG.error(
            f"File name is invalid, use convention [SAMPLE_ID]_R_[1,2].fastq.gz"
        )
        raise click.Abort()

    reference_config = os.path.join(balsamic_cache,
                                    balsamic_version, genome_version,
                                    "reference.json")
    with open(reference_config, 'r') as f:
        reference_dict = json.load(f)["reference"]

    config_collection_dict = PonBalsamicConfigModel(
        QC={
            "quality_trim": quality_trim,
            "adapter_trim": adapter_trim,
            "umi_trim": umi if panel_bed else False,
            "umi_trim_length": umi_trim_length,
        },
        analysis={
            "case_id": case_id,
            "analysis_dir": analysis_dir,
            "analysis_type": "pon",
            "sequencing_type": "targeted" if panel_bed else "wgs"
        },
        reference=reference_dict,
        singularity=os.path.join(balsamic_cache, balsamic_version, "containers"),
        samples=samples,
        panel={
            "capture_kit": panel_bed
        } if panel_bed else None,
    ).dict(by_alias=True, exclude_none=True)
    LOG.info("Config file generated successfully")

    Path.mkdir(Path(config_collection_dict["analysis"]["fastq_path"]),
               parents=True,
               exist_ok=True)
    LOG.info("Directories created successfully")

    create_fastq_symlink(
        casefiles=(normal),
        symlink_dir=Path(config_collection_dict["analysis"]["fastq_path"]),
    )
    LOG.info(f"Symlinks generated successfully")

    config_path = Path(analysis_dir) / case_id / (case_id + "_PON" + ".json")
    with open(config_path, "w+") as fh:
        fh.write(json.dumps(config_collection_dict, indent=4))
    LOG.info(f"Config file saved successfully - {config_path}")

    try:
        generate_graph(config_collection_dict, config_path)
        LOG.info(f"BALSAMIC PON Workflow has been configured successfully!")
    except ValueError as e:
        LOG.error(
            f'BALSAMIC dag graph generation failed - {config_collection_dict["analysis"]["dag"]}',
        )
        raise click.Abort()
