import os
import logging
import glob
import json
import yaml
import click
import shutil

from BALSAMIC.utils.rule import get_result_dir

LOG = logging.getLogger(__name__)


@click.command("scout", short_help="Creates a scout config yaml file.")
@click.option("--sample-config",
              required=True,
              help="Sample config file. Output of balsamic config sample")
@click.pass_context
def scout(context, sample_config):
    '''
    Create a scout config.yaml file
    '''
    LOG.info(f"BALSAMIC started with log level {context.obj['loglevel']}.")
    LOG.info('Adding scout cancer template to delivery directory')

    with open(sample_config, 'r') as fn:
        sample_config = json.load(fn)

    result_dir = get_result_dir(sample_config)
    dst_directory = os.path.join(result_dir, 'delivery_report')
    if not os.path.exists(dst_directory):
        LOG.debug('Creatiing delivery_report directory')
        os.makedirs(dst_directory)

    scout_config_src = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "../../assets/scout_config_template.yaml")
    scout_config_dst = os.path.join(
        dst_directory, sample_config['analysis']['case_id'] + ".scout.yaml")

    LOG.debug('Creating scout config %s', scout_config_dst)
    shutil.copyfile(scout_config_src, scout_config_dst)
    LOG.info('Scout config template is successfully created: %s',
             scout_config_dst)
