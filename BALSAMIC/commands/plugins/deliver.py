import os
import sys
import logging
import glob
import json
import yaml
import click
import snakemake
from collections import defaultdict
from yapf.yapflib.yapf_api import FormatFile

from BALSAMIC.utils.cli import write_json
from BALSAMIC.utils.cli import get_snakefile
from BALSAMIC.utils.cli import CaptureStdout
from BALSAMIC.utils.cli import recursive_default_dict
from BALSAMIC.utils.cli import convert_defaultdict_to_regular_dict
from BALSAMIC.utils.rule import get_result_dir

LOG = logging.getLogger(__name__)


@click.command(
    "deliver",
    short_help=
    "Creates a YAML file with output from variant caller and alignment.")
@click.option("--sample-config",
              required=True,
              help="Sample config file. Output of balsamic config sample")
@click.pass_context
def deliver(context, sample_config):
    '''
    cli for deliver sub-command.
    Writes <case_id>.hk in result_directory.
    '''
    LOG.info(f"BALSAMIC started with log level {context.obj['loglevel']}.")
    LOG.debug("Reading input sample config")
    with open(sample_config, 'r') as fn:
        sample_config_dict = json.load(fn)

    result_dir = get_result_dir(sample_config_dict)
    dst_directory = os.path.join(result_dir, 'delivery_report')
    if not os.path.exists(dst_directory):
        LOG.debug('Creatiing delivery_report directory')
        os.makedirs(dst_directory)
    
    deliver_wildcards = {
        'bam': ['*merged.bam', '*cov.bed'],
        'vcf': ['*.vcf.gz'],
        'vep': ['*.vcf.gz', '*.tsv', '*html', '*balsamic_stat'],
        'cnv': ['*pdf', '*cnr','*cns'],
        'qc': ['multiqc*'],
        'scout': ['*scout.yaml']
    }

    deliver_json = recursive_default_dict()
    for dir_name, file_pattern_list in deliver_wildcards.items():
        deliver_json['files'][dir_name] = list()
        for file_pattern in file_pattern_list:
            list_of_files = glob.glob(
                os.path.join(result_dir, dir_name, file_pattern))
            deliver_json['files'][dir_name].extend(list_of_files)

    yaml_write_directory = os.path.join(result_dir, 'delivery_report')

    os.makedirs(yaml_write_directory, exist_ok=True)

    yaml_file_name = os.path.join(yaml_write_directory,
                                  sample_config_dict['analysis']['case_id'] + ".hk")

    LOG.debug(f"Writing output file {yaml_file_name}")
    deliver_json = convert_defaultdict_to_regular_dict(deliver_json)
    # take the balsamic_stat file
    #balsamic_stat = [s for s in deliver_json['files']['vep'] if s.endswith('balsamic_stat')][0]
    #with open(balsamic_stat, 'r') as f:
    #    balsamic_stat = yaml.load(f, Loader=yaml.SafeLoader) 

    #deliver_json['meta'] = balsamic_stat

#    with open(yaml_file_name, 'w') as f:
#        yaml.dump(deliver_json, f, default_flow_style=False)

    #yaml.dump(deliver_json, sys.stdout, default_flow_style=False)
    analysis_type = sample_config_dict['analysis']['analysis_type']
    sequencing_type = sample_config_dict['analysis']['sequencing_type']
    snakefile = get_snakefile(analysis_type, sequencing_type)
    delivery_file_name = os.path.join(yaml_write_directory,
                                  sample_config_dict['analysis']['case_id'] + ".hk")
    delivery_file_raw = os.path.join(yaml_write_directory,
                                  sample_config_dict['analysis']['case_id'] + "_delivery_raw.hk")
    with open(delivery_file_raw, 'r') as fn:
        delivery_file_raw_dict = json.load(fn)

    snakemake.snakemake(snakefile=snakefile, config={'delivery':'True'}, dryrun=True, summary=True, configfile=sample_config, quiet=True)
    with CaptureStdout() as summary:
        snakemake.snakemake(snakefile=snakefile, config={'delivery':'True'}, dryrun=True, summary=True, configfile=sample_config, quiet=True)
    summary = [i.split('\t') for i in summary]
    summary_dict = [dict(zip(summary[0], value)) for value in summary[1:]]

    output_files_merged = defaultdict(dict)
    for interm_list in (summary_dict, delivery_file_raw_dict):
        for item in interm_list:
            output_files_merged[item['output_file']].update(item)

    delivery_json = dict()
    delivery_json['files'] = list()
    for k,v in output_files_merged.items():
        if 'date' in output_files_merged[k]:
            delivery_json['files'].append(output_files_merged[k])
#    output_files_merged = output_files_merged.values()
    write_json(delivery_json, delivery_file_name)
    #yaml.dump([dict(zip(summary[0], value)) for value in summary[1:]], sys.stdout, default_flow_style=False)













