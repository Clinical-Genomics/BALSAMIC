import os
import re
import json
import shutil
import logging
import click
import snakemake
import graphviz
import sys
import yaml
import BALSAMIC

from pathlib import Path
from datetime import datetime
from yapf.yapflib.yapf_api import FormatFile
from BALSAMIC.utils.cli import CaptureStdout
from BALSAMIC.utils.errors import ExtensionNotSupportedError
from BALSAMIC.config.constats import *

LOG = logging.getLogger(__name__)

CASE_COMMAND = click.command(
    "case",
    short_help="Create a sample config file from input sample data")

CASE_ID_OPTION = click.Option(
    "--case-id",
    required=True,
    help="Sample id that is used for reporting, \
    naming the analysis jobs, and analysis path"
    )

UMI_OPTION = click.Option(
    "--umi/--no-umi",
    default=True,
    show_default=True,
    is_flag=True,
    help="UMI processing steps for samples with UMI tags"
    )
UMI_TRIM_LENGTH_OPTION = click.Option(
    "--umi-trim-length",
    default=5,
    show_default=True,
    type=int,
    help="Trim N bases from reads in fastq"
    )
QUALITY_TRIM_OPTION = click.Option(
    "--quality-trim/--no-quality-trim",
    default=True,
    show_default=True,
    is_flag=True,
    help="Trim low quality reads in fastq"
    )
ADAPTER_TRIM_OPTION = click.Option(
    "--adapter-trim/--no-adapter-trim",
    default=False,
    show_default=True,
    is_flag=True,
    help="Trim adapters from reads in fastq"
    )
REFERENCE_CONFIG_OPTION = click.Option(
    "-r", "--reference-config",
    required=True,
    type=click.Path(),
    exists=True,
    resolve_path=True,
    help="Reference config file."
    )
PANEL_BED_OPTION = click.Option(
    "-p","--panel-bed",
    type=click.Path(),
    required=False,
    exists=True,
    resolve_path=True,
    help="Panel bed file for variant calling."
    )
SINGULARITY_OPTION = click.Option(
    "--singularity",
    type=click.Path(),
    required=True,
    resolve_path=True,
    exists=True,
    help="Download singularity image for BALSAMIC"
    )
FASTQ_PREFIX_OPTION = click.Option(
    "--fastq-prefix",
    required=False,
    default="",
    help="Prefix to fastq file. The string that comes after readprefix"
    )
ANALYSIS_DIR_OPTION = click.Option(
    "--analysis-dir",
    type=click.Path(),
    default=".",
    resolve_path=True,
    exists=True,
    help="Root analysis path to store analysis logs and results. \
                                     The final path will be analysis-dir/sample-id"
    )
TUMOR_OPTION = click.Option(
    "-t","--tumor",
    required=True,
    multiple=True,
    exists=True,
    resolve_path=True,
    help="Fastq files for tumor sample."
    )
NORMAL_OPTION = click.Option(
    "-n","--normal",
    required=False,
    multiple=True,
    exists=True,
    resolve_path=True,
    help="Fastq files for normal sample."
    )
FORMAT_OPTION = click.Option(
    "-f", "--format", 
    required=False, 
    type=click.Choice(["fastq", "vcf", "bam"]), 
    default="fastq"
    )


def validate_fastq_pattern(sample):
    fq_pattern = re.compile(r"R_[12]"+".fastq.gz$")
    sample_basename = os.path.basename(sample)
    try:
        file_str = sample_basename[0:(
            fq_pattern.search(sample_basename).span()[0] + 1)]
        return file_str

    except AttributeError:
        LOG.error(
            f"File name is invalid, fastq file should be sample_R_1.fastq.gz")


def get_panel_chrom(panel_bed) -> list:
    lines = [line.rstrip('\n') for line in open(panel_bed, 'r')]
    return list(set([s.split('\t')[0] for s in lines]))
    
def get_bioinfo_tools_list(conda_env_path) -> dict:
    bioinfo_tools = {}
    for yaml_file in os.listdir(conda_env_path):
        if yaml_file.endswith(".yaml"):
            with open(conda_env_path + yaml_file, "r") as f:
                packages = yaml.safe_load(f)["dependencies"]
                for p in packages:
                    try:
                        name, version = p.split("=")
                    except ValueError:
                        name, version = p, ""
                    finally:
                        if name in core_packages:
                            bioinfo_tools[name] = version
    return bioinfo_tools


def get_sample_dict(tumor, normal):
    samples = {}
    if normal:
        for sample in normal:
            key, val = get_sample_names(sample, "normal")
            samples[key] = val

    for sample in tumor:
        key, val = get_sample_names(sample, "tumor")
        samples[key] = val

    return samples
    

def get_sample_names(file, sample_type):
    file_str = validate_fastq_pattern(file)
    return file_str, {
        "file_prefix": file_str,
        "type": sample_type,
        "readpair_suffix": ["1", "2"]}


@CASE_COMMAND
@CASE_ID_OPTION
@UMI_OPTION
@UMI_TRIM_LENGTH_OPTION
@QUALITY_TRIM_OPTION
@ADAPTER_TRIM_OPTION
@REFERENCE_CONFIG_OPTION
@PANEL_BED_OPTION
@SINGULARITY_OPTION
@ANALYSIS_DIR_OPTION
@TUMOR_OPTION
@NORMAL_OPTION
@FORMAT_OPTION
@click.pass_context
def case_config(context, case_id, umi, umi_trim_length, adapter_trim, quality_trim, reference_config, panel_bed, singularity, analysis_dir, tumor, normal, format):

    reference_dict = json.load(open(reference_config))["reference"]
    bioinfo_tools = get_bioinfo_tools_list(CONDA_ENV_PATH)
    samples = get_sample_dict(tumor, normal)

    config_collection = BalsamicConfigModel(QC={
                                                "quality_trim" : quality_trim,
                                                "adapter_trim" : adapter_trim,
                                                "umi_trim" : umi,
                                                "umi_trim_length" : umi_trim_length,
                                                },
                                            analysis={
                                                "case_id": case_id,
                                                "analysis_dir": analysis_dir,
                                                "analysis_type" : "paired" if normal else "single",
                                                "sequencing_type" : "targeted" if panel_bed else "wgs",
                                                },
                                            bioinfo_tools=bioinfo_tools,
                                            samples=samples,
                                            reference=reference_dict,
                                            panel={
                                                "panel" : panel_bed if panel_bed else None, 
                                                "chrom" : get_panel_chrom(panel_bed) if panel_bed else None, 
                                                },

                                            singularity=singularity,
                                            )
