"""
Entry cli for balsamic
"""
import click

# Subcommands
from BALSAMIC.commands.install import install as install_command
from BALSAMIC.commands.run import run_analysis as analysis_command
from BALSAMIC.commands.config import config as config_command
from BALSAMIC.commands.report import report as report_command

# CLI commands and decorators
from BALSAMIC.tools.cli_utils import add_doc as doc

# Get version
from BALSAMIC import __version__


@click.group()
@click.version_option(version=__version__)
@click.pass_context
@doc("""BALSAMIC {version}: Bioinformatic Analysis pipeLine for
        SomAtic MutatIons in Cancer""".format(version=__version__))
def cli():
    "BALSAMIC"


cli.add_command(install_command)
cli.add_command(analysis_command)
cli.add_command(config_command)
cli.add_command(report_command)
