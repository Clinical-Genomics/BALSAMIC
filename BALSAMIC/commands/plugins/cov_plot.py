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
@click.option("--sample-config",
              required=True,
              help="Sample config file. Output of balsamic config sample")
@click.pass_context
def target_cov_plot(context, sample_config):
    '''
    cli for coverage plot sub-command.
    Creates coverage plots in result_directory.
    '''

    LOG.debug("Reading input sample config")
    with open(sample_config, 'r') as fn:
        sample_config = json.load(fn)


    result_dir = get_result_dir(sample_config)
    # List of bam files:
    cov_files = glob.glob(os.path.join(result_dir, "bam", "*.cov.bed"))

    for f in cov_files:
        LOG.info("Plotting coverage for %s", f)
#        plot_file = plot_cov.hist(cov_file=f)
        LOG.info("Coverage plot file: %s", plot_file)
# TODO: patten for an input region
# pattern = re.compile("^([1-9]{1,2}):([0-9]+)(\-)([0-9]+)") 

    LOG.debug("Plotting coverage")
