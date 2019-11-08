import os
import sys
import logging
import glob
import json
import yaml
import click

from BALSAMIC.utils.rule import get_result_dir

LOG = logging.getLogger(__name__)


@click.command(
    "housekeeper",
    short_help=
    "Creates a YAML file with output from variant caller and alignment.")
@click.option("--sample-config",
              required=True,
              help="Sample config file. Output of balsamic config sample")
@click.pass_context
def housekeeper(context, sample_config):
    '''
    cli for housekeeper sub-command.
    Writes <case_id>.hk in result_directory.
    '''
    LOG.info(f"BALSAMIC started with log level {context.obj['loglevel']}.")
    LOG.debug("Reading input sample config")
    with open(sample_config, 'r') as fn:
        sample_config = json.load(fn)

    result_dir = get_result_dir(sample_config)
    dst_directory = os.path.join(result_dir, 'delivery_report')
    if not os.path.exists(dst_directory):
        LOG.debug('Creatiing delivery_report directory')
        os.makedirs(dst_directory)
    
    housekeeper_wildcards = {
        'bam': ['*merged.bam', '*cov.bed'],
        'vcf': ['*.vcf.gz'],
        'vep': ['*.vcf.gz'],
        'cnv': ['*pdf', '*cnr','*cns'],
        'qc': ['multiqc*'],
        'scout': ['*scout.yaml']
    }

    housekeeper_json = dict()
    for dir_name, file_pattern_list in housekeeper_wildcards.items():
        housekeeper_json[dir_name] = list()
        for file_pattern in file_pattern_list:
            list_of_files = glob.glob(
                os.path.join(result_dir, dir_name, file_pattern))
            housekeeper_json[dir_name].extend(list_of_files)

    yaml_write_directory = os.path.join(result_dir, 'delivery_report')

    os.makedirs(yaml_write_directory, exist_ok=True)

    yaml_file_name = os.path.join(yaml_write_directory,
                                  sample_config['analysis']['case_id'] + ".hk")

    LOG.debug(f"Writing output file {yaml_file_name}")
    with open(yaml_file_name, 'w') as f:
        yaml.dump(housekeeper_json, f, default_flow_style=False)

    yaml.dump(housekeeper_json, sys.stdout, default_flow_style=False)
