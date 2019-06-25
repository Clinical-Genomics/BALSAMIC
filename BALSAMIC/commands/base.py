"""
Entry cli for balsamic
"""
import click

# Subcommands
from BALSAMIC.commands.initiate.initiate import initiate as initiate_command
from BALSAMIC.commands.install.install import install as install_command
from BALSAMIC.commands.run.run_analysis import run_analysis as analysis_command
from BALSAMIC.commands.config.base import config as config_command

# CLI commands and decorators
from BALSAMIC.utils.cli import add_doc as doc

# Get version
from BALSAMIC import __version__

@click.group()
@click.version_option(version=__version__)
@click.pass_context

@doc("""BALSAMIC {version}: Bioinformatic Analysis pipeLine for
        SomAtic MutatIons in Cancer""".format(version=__version__))
def cli(context):
    "BALSAMIC"

cli.add_command(initiate_command)
cli.add_command(install_command)
cli.add_command(analysis_command)
cli.add_command(config_command)
