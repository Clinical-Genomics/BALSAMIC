import os
import logging
import glob
import json
import yaml
import click

from BALSAMIC.utils.rule import get_result_dir

LOG = logging.getLogger(__name__)


@click.command(
    "scout",
    short_help=
    "Creates a scout config yaml file.")
@click.option("--sample-config",
              required=True,
              help="Sample config file. Output of balsamic config sample")
@click.pass_context
def scout(context, sample_config):
    '''
    Create a scout config.yaml file
    '''
    
    LOG.info('Placeholder for possible scout config generation')
