"""Balsamic configuration file generation base commands."""
import click

from BALSAMIC.commands.config.case import case_config
from BALSAMIC.commands.config.pon import pon_config


@click.group()
@click.pass_context
def config():
    """Create config files required for running the pipeline."""
    pass


config.add_command(case_config)
config.add_command(pon_config)
