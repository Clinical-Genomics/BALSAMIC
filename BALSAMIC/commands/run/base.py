"""Balsamic run CLI."""
import click

from BALSAMIC.commands.run.analysis import analysis as run_analysis_command


@click.group()
@click.pass_context
def run(context: click.Context):
    """Run Balsamic analysis on a provided configuration file."""
    pass


run.add_command(run_analysis_command)
