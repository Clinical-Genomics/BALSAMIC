#!/usr/bin/env python
import os
import logging
import glob
import json
import click

from BALSAMIC.utils.rule import get_result_dir

LOG = logging.getLogger(__name__)

@click.command(
    "delivery",
    short_help=
    "Creates a YAML file with output from variant caller and alignment.")
@click.option(
    "--sample-config",
    required=True,
    help="Sample config file. Output of balsamic config sample")
@click.option(
    "--output-dir",
    required=True,
    help="Output directory for writing the json file.")
@click.pass_context
def delivery(context, sample_config, output_dir):
    '''
    cli for delivery sub-command
    '''
    
    LOG.debug("Reading input sample config")
    with open(sample_config, 'r') as fn:
        sample_config = json.load(fn)

    delivery_wildcards = {'bam':['*merged.bam'], 'vcf':['*.vcf.gz'], 'vep':['*.vcf.gz', '*.txt']}

    delivery_json = dict()
    for dir_name, file_pattern_list in delivery_wildcards.items():
        delivery_json[dir_name] = list()
        result_dir = os.path.join(get_result_dir(sample_config), dir_name)
        for file_pattern in file_pattern_list:
            list_of_files = glob.glob(os.path.join(result_dir, file_pattern))
            delivery_json[dir_name].extend(list_of_files)

    print(json.dumps(delivery_json, indent=4))
