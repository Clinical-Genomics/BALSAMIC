import re
import os
import logging
import glob
import json
import yaml
import click

from BALSAMIC.utils.rule import get_result_dir
#from BALSAMIC.utils import plot_cov

LOG = logging.getLogger(__name__)


@click.command(
    "target-cov-plot",
    short_help=
    "Plots coverage for target regions.")
@click.option("--input-cov",
              required=True,
              help="Sample config file. Output of balsamic config sample")
@click.pass_context
def target_cov_plot(context, input_cov):
    '''
    cli for coverage plot sub-command.
    Creates coverage plots in result_directory.
    '''
    
    

