#!/usr/bin/env python
import os
import subprocess
import hashlib
import yaml
import click
import logging
import sys
import json
from yapf.yapflib.yapf_api import FormatFile
from datetime import datetime

def merge_json(*args):
    '''
    Take a list of json files and merges them together

    Input: list of json file
    Output: dictionary of merged json
    '''

    json_out = dict()
    for json_file in args:
        try:
            if isinstance(json_file, dict):
                json_out = {**json_out, **json_file}
            else:
                with open(json_file) as fn:
                    json_out = {**json_out, **json.load(fn)}
        except OSError:
            print("File not found")

    return json_out


def write_json(json_out, output_config):

    try:
        with open(output_config, "w") as fn:
            json.dump(json_out, fn)
    except OSError:
        print("Write failed")


def get_config(config_name):
    """
    Return a string path for config file.
    """

    try:
        config_file = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), config_name + '.json')
    except OSError:
        print("Couldn't locate config file" + config_name + ".json")

    return config_file


def set_panel_bed(json_out, panel_bed):
    """
    Set panel path in config file
    """
    try:
        json_out['path']['panel'] = os.path.split(
            os.path.abspath(panel_bed))[0] + "/"
        json_out['bed']['variant_panel'] = os.path.split(
            os.path.abspath(panel_bed))[1]
    except OSError:
        print("Couldn't locate bed file" + panel_bed)

    return json_out


def check_exist(path):
    """
    Checks if fastq file readable and accessable.
    """

    try:
        f = open(path, 'r')
        f.close()
    except (IOError, FileNotFoundError) as e:
        print("File not found or unreadable.", path)
        raise e

    return True


@click.command(
    "sample", short_help="Create a sample config file from input sample data")
@click.option(
    '-a',
    '--analysis-type',
    required=False,
    default='paired',
    show_default=True,
    type=click.Choice(['paired', 'single']),
    help=
    "Analysis config file for paired (tumor vs normal) or single (tumor-only) mode."
)
@click.option(
    '-i',
    '--install-config',
    required=False,
    default=get_config("install"),
    show_default=True,
    type=click.Path(),
    help="Installation config file.")
@click.option(
    '-r',
    '--reference-config',
    required=False,
    default=get_config("reference"),
    show_default=True,
    type=click.Path(),
    help='Reference config file.')
@click.option(
    '-p',
    '--panel-bed',
    required=True,
    type=click.Path(),
    help='Panel bed file for variant calling.')
@click.option(
    '-s',
    '--sample-config',
    type=click.Path(),
    help='Input sample config file.')
@click.option(
    '-o',
    '--output-config',
    required=True,
    type=click.Path(),
    help='Output a json config file ready to be imported for run-analysis')
@click.option(
    '-t',
    '--tumor',
    required=True,
    help=
    'Fastq files for tumor sample. Example: if files are tumor_fqreads_1.fastq.gz tumor_fqreads_2.fastq.gz, the input should be --tumor tumor_fqreads'
)
@click.option(
    '-n',
    '--normal',
    help=
    'Fastq files for normal sample. Example: if files are normal_fqreads_1.fastq.gz normal_fqreads_2.fastq.gz, the input should be --normal normal_fqreads'
)
@click.option(
    '--sample-id',
    help=
    'Sample id that is used for reporting, naming the analysis jobs, and analysis path'
)
@click.option(
    '--analysis-dir',
    type=click.Path(),
    help=
    'Root analysis path to store analysis logs and results. The final path will be analysis-dir/sample-id'
)
@click.option(
    '--fastq-path',
    type=click.Path(),
    help=
    'Path for fastq files. All fastq files should be within same path and that path has to exist.'
)
@click.option(
    '--check-files/--no-check-files',
    default=True,
    help='Check if fastq input files exist')
@click.pass_context
def sample(context, analysis_type, install_config, sample_config,
           reference_config, panel_bed, output_config, normal, tumor,
           sample_id, analysis_dir, fastq_path, check_files):
    """
    Prepares a config file for balsamic run_analysis. For now it is just treating json as dictionary and merging them as
it is. So this is just a placeholder for future.

    """
    analysis_config = get_config("analysis_" + analysis_type)

    click.echo("Reading analysis config file %s" % analysis_config)
    click.echo("Reading reference config file %s" % reference_config)

    read_prefix = ["1", "2"]

    if sample_config:
        click.echo("Reading sample config file %s" % sample_config)
        sample_config = os.path.abspath(sample_config)
    else:
        sample_config = get_config("sample")
        click.echo("Reading sample config file %s" % sample_config)
        with open(sample_config) as j:
            sample_config = json.load(j)
        sample_config["analysis"]["sample_id"] = sample_id
        sample_config["analysis"][
            "config_creation_date"] = datetime.now().strftime("%Y-%m-%d %H:%M")
        sample_config["analysis"]["analysis_dir"] = analysis_dir + "/"
        sample_config["analysis"]["fastq_path"] = fastq_path + "/"
        sample_config["analysis"]["analysis_type"] = analysis_type
        sample_config["samples"] = {}

        sample_config["samples"][tumor] = {
            "file_prefix": tumor,
            "type": "tumor",
            "readpair_suffix": read_prefix
        }
        if check_files:
            fq_list = [
                os.path.join(fastq_path, tumor + "_" + r + ".fastq.gz")
                for r in read_prefix
            ]
            for f in fq_list:
                check_exist(f)

        if normal:
            sample_config["samples"][normal] = {
                "file_prefix": normal,
                "type": "normal",
                "readpair_suffix": read_prefix
            }
            if check_files:
                fq_list = [
                    os.path.join(fastq_path, tumor + "_" + r + ".fastq.gz")
                    for r in read_prefix
                ]
                for f in fq_list:
                    check_exist(f)

    json_out = merge_json(analysis_config, sample_config, reference_config,
                          install_config)

    output_config = os.path.join(
        os.path.abspath(json_out["analysis"]["analysis_dir"]),
        json_out["analysis"]["sample_id"], output_config)
    click.echo("Writing output config file %s" % output_config)

    if panel_bed:
        json_out = set_panel_bed(json_out, panel_bed)

    os.makedirs(os.path.dirname(os.path.abspath(output_config)), exist_ok=True)
    write_json(json_out, output_config)
    FormatFile(output_config, in_place=True)
