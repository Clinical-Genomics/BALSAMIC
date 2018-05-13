#!/usr/bin/env python
import os
import subprocess
import hashlib
import yaml
import click
import logging
import sys
import json


def merge_json(*args):
    '''
    Take a list of json files and merges them together

    Input: list of json file
    Output: dictionary of merged json
    '''

    json_out = dict()
    for json_file in args:
        print(json_file)
        try:
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
        json_out['path']['panel'] = os.path.split(os.path.abspath(panel_bed))[0] + "/"
        json_out['bed']['variant_panel'] = os.path.split(os.path.abspath(panel_bed))[1]
    except OSError:
        print("Couldn't locate bed file" + panel_bed)

    return json_out

@click.command(
    "create_config",
    short_help="Create a sample config file from input sample data")
@click.option(
    '-a',
    '--analysis-config',
    required=False,
    default=get_config("analysis"),
    show_default=True,
    type=click.Path(),
    help="Analysis config file.")
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
    required=True,
    type=click.Path(),
    help='Input sample config file.')
@click.option(
    '-o',
    '--output-config',
    required=True,
    type=click.Path(),
    help='Output a json config file ready to be imported for run-analysis')
@click.pass_context
def create_config(context, analysis_config, install_config, sample_config,
                  reference_config, panel_bed, output_config):
    """
    Prepares a config file for balsamic run_analysis. For now it is just treating json as dictionary and merging them as
it is. So this is just a placeholder for future.

    """

    click.echo("Reading analysis config file %s" % analysis_config)
    click.echo("Reading sample config file %s" % sample_config)
    click.echo("Reading reference config file %s" % reference_config)
    click.echo("Writing output config file %s" % output_config)

    json_out = merge_json(analysis_config, sample_config, reference_config,
                          install_config)
    
    if panel_bed:
        json_out = set_panel_bed(json_out, panel_bed)

    write_json(json_out, output_config)
