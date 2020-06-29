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
from BALSAMIC.config.constats import QCModel, VCFModel, AnalysisModel, SampleInstanceModel, BioinfoToolsModel, PanelModel, BalsamicConfigModel, CONDA_ENV_PATH, CONDA_ENV_YAML, RULE_DIRECTORY

LOG = logging.getLogger(__name__)


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
    for yaml_file in Path(conda_env_path).rglob('*.yaml'):
        with open(yaml_file, "r") as f:
            packages = yaml.safe_load(f)["dependencies"]
            for p in packages:
                try:
                    name, version = p.split("=")
                except ValueError:
                    name, version = p, ""
                finally:
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



@click.command("case",
    short_help="Create a sample config file from input sample data")
@click.option(
    "--case-id",
    required=True,
    help="Sample id that is used for reporting, \
    naming the analysis jobs, and analysis path"
    )
@click.option(
    "--umi/--no-umi",
    default=True,
    show_default=True,
    is_flag=True,
    help="UMI processing steps for samples with UMI tags"
    )
@click.option(
    "--umi-trim-length",
    default=5,
    show_default=True,
    type=int,
    help="Trim N bases from reads in fastq"
    )
@click.option(
    "--quality-trim/--no-quality-trim",
    default=True,
    show_default=True,
    is_flag=True,
    help="Trim low quality reads in fastq"
    )
@click.option(
    "--adapter-trim/--no-adapter-trim",
    default=False,
    show_default=True,
    is_flag=True,
    help="Trim adapters from reads in fastq"
    )
@click.option(
    "-r", "--reference-config",
    required=True,
    type=click.Path(exists=True, resolve_path=True),
    help="Reference config file."
    )
@click.option(
    "-p","--panel-bed",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    help="Panel bed file for variant calling."
    )
@click.option(
    "--singularity",
    type=click.Path(exists=True, resolve_path=True),
    required=True,
    help="Download singularity image for BALSAMIC"
    )
@click.option(
    "--fastq-prefix",
    required=False,
    default="",
    help="Prefix to fastq file. The string that comes after readprefix"
    )
@click.option(
    "--analysis-dir",
    type=click.Path(exists=True, resolve_path=True),
    default=".",
    help="Root analysis path to store analysis logs and results. \
                                     The final path will be analysis-dir/sample-id"
    )
@click.option(
    "-t","--tumor",
    type=click.Path(exists=True, resolve_path=True),
    required=True,
    multiple=True,
    help="Fastq files for tumor sample."
    )
@click.option(
    "-n","--normal",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    multiple=True,
    help="Fastq files for normal sample."
    )
@click.option(
    "-f", "--format", 
    required=False, 
    type=click.Choice(["fastq", "vcf", "bam"]), 
    default="fastq",
    show_default=True,
    )
@click.pass_context
def case_config(context, case_id, umi, umi_trim_length, adapter_trim, quality_trim, fastq_prefix, reference_config, panel_bed, singularity, analysis_dir, tumor, normal, format):

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
                                            panel={
                                                "capture_kit" : panel_bed, 
                                                "chrom" : get_panel_chrom(panel_bed), 
                                                } if panel_bed else None,
                                            bioinfo_tools=bioinfo_tools,
                                            reference=reference_dict,
                                            singularity=singularity,
                                            samples=samples,
                                            vcf={},
                                            )
    print(json.dumps(config_collection.dict(), indent=4))
    #print(config_collection.dict())

    #Make folders
    #Create Symlinks
    #Save json
    #Create DAG
    

if __name__ == "__main__":
    case_config()