#!/usr/bin/env python
import os
import subprocess
import re
import json
import copy
import glob
import click
import snakemake
import graphviz
from datetime import datetime
from pathlib import Path
from yapf.yapflib.yapf_api import FormatFile

from BALSAMIC.utils.cli import get_package_split
from BALSAMIC.utils.cli import write_json
from BALSAMIC.utils.cli import get_config
from BALSAMIC.utils.cli import get_ref_path
from BALSAMIC.utils.cli import get_snakefile
from BALSAMIC.utils.cli import CaptureStdout
from BALSAMIC.utils.rule import get_chrom
from BALSAMIC import __version__ as bv


def merge_json(*args):
    """
    Take a list of json files and merges them together

    Input: list of json file
    Output: dictionary of merged json
    """

    json_out = dict()
    for json_file in args:
        try:
            if isinstance(json_file, dict):
                json_out = {**json_out, **json_file}
            else:
                with open(json_file) as fn:
                    json_out = {**json_out, **json.load(fn)}
        except OSError as error:
            raise error

    return json_out


def set_panel_bed(json_out, panel_bed):
    """
    Set panel path in config file
    """
    try:
        json_out["panel"] = {"capture_kit": os.path.abspath(panel_bed)}
        json_out["panel"]["chrom"] = get_chrom(panel_bed)

    except OSError as error:
        raise error

    return json_out


def check_exist(path):
    """
    Checks if fastq file readable and accessable.
    """

    try:
        f = open(path, "r")
        f.close()
    except (IOError, FileNotFoundError) as error:
        raise error

    return True


def get_analysis_type(normal, umi):
    """ return analysis type """
    if umi:
        return "paired_umi" if normal else "single_umi"

    return "paired" if normal else "single"


def get_output_config(config, sample_id):
    """ return output config json file"""
    if not config:
        return sample_id + "_" + datetime.now().strftime("%Y%m%d") + ".json"
    else:
        return config


def get_sample_config(sample_config, sample_id, analysis_dir, analysis_type):
    """
    creating sample config to run the analysis
    """
    with open(sample_config) as sample_json:
        sample_config = json.load(sample_json)

    sample_config["analysis"]["sample_id"] = sample_id
    sample_config["analysis"]["config_creation_date"] = datetime.now(
    ).strftime("%Y-%m-%d %H:%M")
    sample_config["analysis"]["analysis_dir"] = analysis_dir + "/"
    sample_config["analysis"]["log"] = os.path.join(analysis_dir, sample_id,
                                                    'logs/')
    sample_config["analysis"]["script"] = os.path.join(analysis_dir, sample_id,
                                                       'scripts/')
    sample_config["analysis"]["result"] = os.path.join(analysis_dir, sample_id,
                                                       'analysis')
    sample_config["analysis"]["analysis_type"] = analysis_type
    sample_config["samples"] = {}

    return sample_config


def link_fastq(src_files, des_path):
    """
    Creating fastq symlinks in given destination
    """
    for src_file in src_files:
        basename = os.path.basename(src_file)
        des_file = os.path.join(des_path, basename)
        try:
            os.symlink(src_file, des_file)
        except FileExistsError:
            print(
                f"Desitination file {des_file} exists. No symbolic link was created."
            )


def get_fastq_path(file, fq_pattern):
    # check the fastq if exists
    file = os.path.abspath(file)
    if Path(file).exists():
        file_basename = os.path.basename(file)
        try:
            # extracting file prefix
            file_str = file_basename[0:(fq_pattern.search(file_basename).span()[0] + 1)]
        except AttributeError as error:
            raise error
    else:
        raise Exception(FileNotFoundError)

    return file_str, os.path.split(file)[0]


def configure_fastq(fq_path, sample, fastq_prefix):
    """
    Configure the fastq files for analysis
    """
    fq_pattern = re.compile(r"R_[12]" + fastq_prefix + ".fastq.gz$")
    paths = list()

    # get a list of fq files
    sample_str, sample_path = get_fastq_path(sample, fq_pattern)
    paths.append(sample_path)

#    if normal:
#        normal_str, normal_path = get_fastq_path(normal, fq_pattern)
#        paths.append(normal_path)

    fq_files = set()
    for path in paths:
        for file in os.listdir(path):
            if fq_pattern.search(file):
                fq_files.add(os.path.join(path, file))

    # create symlink
    link_fastq(fq_files, fq_path)

    # return file prefix
    return sample_str 


@click.command("sample",
               short_help="Create a sample config file from input sample data")
@click.option('--umi',
              is_flag=True,
              help="UMI processing steps for samples with umi tags")
@click.option("-i",
              "--install-config",
              required=False,
              default=get_config("install"),
              show_default=True,
              type=click.Path(),
              help="Installation config file.")
@click.option("-r",
              "--reference-config",
              required=True,
              show_default=True,
              type=click.Path(),
              help="Reference config file.")
@click.option("-p",
              "--panel-bed",
              required=True,
              type=click.Path(),
              help="Panel bed file for variant calling.")
@click.option("-s",
              "--sample-config",
              type=click.Path(),
              help="Input sample config file.")
@click.option(
    "-o",
    "--output-config",
    required=False,
    help="Output a json config filename ready to be imported for run-analysis")
@click.option(
    "-t",
    "--tumor",
    required=True,
    help="Fastq files for tumor sample. \
              Example: if files are tumor_fqreads_1.fastq.gz tumor_fqreads_2.fastq.gz, \
              the input should be --tumor tumor_fqreads",
)
@click.option("-n",
              "--normal",
              help="Fastq files for normal sample. \
              Example: if files are normal_fqreads_1.fastq.gz normal_fqreads_2.fastq.gz, \
              the input should be --normal normal_fqreads")
@click.option("--sample-id",
              required=True,
              help="Sample id that is used for reporting, \
              naming the analysis jobs, and analysis path")
@click.option("--fastq-prefix",
              required=False,
              default="",
              help="Prefix to fastq file. \
              The string that comes after readprefix")
@click.option("--analysis-dir",
              type=click.Path(),
              help="Root analysis path to store \
              analysis logs and results. The final path will be analysis-dir/sample-id"
              )
@click.option("--overwrite-config/--no-overwrite-config",
              default=True,
              help="Overwrite output config file")
@click.option("--create-dir/--no-create-dir",
              default=True,
              help="Create analysis directiry.")
@click.pass_context
def sample(context, umi, install_config, sample_config, reference_config,
           panel_bed, output_config, normal, tumor, sample_id, analysis_dir,
           overwrite_config, create_dir, fastq_prefix):
    """
    Prepares a config file for balsamic run_analysis. For now it is just treating json as
    dictionary and merging them as it is. So this is just a placeholder for future.
    """

    analysis_type = get_analysis_type(normal, umi)
    output_config = get_output_config(output_config, sample_id)
    analysis_config = get_config("analysis_" + analysis_type)

    click.echo("Reading analysis config file %s" % analysis_config)
    click.echo("Reading reference config file %s" % reference_config)

    reference_json = get_ref_path(reference_config)

    read_prefix = ["1", "2"]

    if sample_config:
        sample_config_path = os.path.abspath(sample_config)
    else:
        sample_config_path = get_config("sample")

    click.echo("Reading sample config file %s" % sample_config_path)

    analysis_dir = os.path.abspath(analysis_dir)
    sample_config = get_sample_config(sample_config_path, sample_id,
                                      analysis_dir, analysis_type)

    output_dir = os.path.join(analysis_dir, sample_id)

    if create_dir:
        os.makedirs(output_dir, exist_ok=True)

    # create dir for fastq symlink creation
    fq_path = os.path.join(output_dir, 'analysis', 'fastq')
    os.makedirs(fq_path, exist_ok=True)

    tumor = configure_fastq(fq_path, tumor, fastq_prefix)
    if normal:
        normal = configure_fastq(fq_path, normal, fastq_prefix)

    sample_config["samples"][tumor] = {
        "file_prefix": tumor,
        "type": "tumor",
        "readpair_suffix": read_prefix,
    }

    if normal:
        sample_config["samples"][normal] = {
            "file_prefix": normal,
            "type": "normal",
            "readpair_suffix": read_prefix,
        }

    sample_config["analysis"]["fastq_path"] = fq_path + "/"
    sample_config["analysis"]["BALSAMIC_version"] = bv

    conda_env = glob.glob(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "../..",
                     "conda/*.yaml"))

    bioinfo_config = dict()
    bioinfo_config["bioinfo_tools"] = get_package_split(conda_env)

    output_config = os.path.join(output_dir, output_config)
    click.echo(
        "Writing output config file %s" % os.path.abspath(output_config))

    json_out = merge_json(analysis_config, sample_config, reference_json,
                          install_config, bioinfo_config)

    dag_image = os.path.join(output_dir,
                             output_config + '_BALSAMIC_' + bv + '_graph')

    json_out["analysis"]["dag"] = dag_image + ".pdf"

    if panel_bed:
        json_out = set_panel_bed(json_out, panel_bed)

    if overwrite_config:
        write_json(json_out, output_config)

    FormatFile(output_config, in_place=True)

    with CaptureStdout() as graph_dot:
        snakemake.snakemake(snakefile=get_snakefile(analysis_type),
                            dryrun=True,
                            configfile=output_config,
                            printrulegraph=True)

    graph_title = "_".join(['BALSAMIC', bv, json_out["analysis"]["sample_id"]])
    graph_dot = "".join(graph_dot).replace(
        'snakemake_dag {',
        'BALSAMIC { label="' + graph_title + '";labelloc="t";')
    graph_obj = graphviz.Source(graph_dot,
                                filename=dag_image,
                                format="pdf",
                                engine="dot")
    #    graph_obj.attr('graph',label='BALSAMIC')
    #    graph_obj.graph_attr['label'] = "_".join(['BALSAMIC',bv,json_out["analysis"]["sample_id"]])
    graph_obj.render()
