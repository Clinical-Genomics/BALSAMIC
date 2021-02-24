import logging
import click

LOG = logging.getLogger(__name__)


@click.command("target-cov-plot",
               short_help="Plots coverage for target regions.")
@click.pass_context
def target_cov_plot(context):
    '''
    cli for coverage plot sub-command.
    Creates coverage plots in result_directory.
    '''
