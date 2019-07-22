"""
Entry cli for balsamic
"""
import logging
import click
import coloredlogs

# Subcommands
from BALSAMIC.commands.install.install import install as install_command
from BALSAMIC.commands.run.base import run as run_command
from BALSAMIC.commands.config.base import config as config_command

# CLI commands and decorators
from BALSAMIC.utils.cli import add_doc as doc

# Get version
from BALSAMIC import __version__

LOG = logging.getLogger(__name__)


@click.group()
@click.version_option(version=__version__)
@click.pass_context
@doc("""BALSAMIC {version}: Bioinformatic Analysis pipeLine for
        SomAtic MutatIons in Cancer""".format(version=__version__))
def cli(context):
    "BALSAMIC"
    coloredlogs.install(level='INFO', fmt='%(asctime)s %(hostname)s %(name)s %(levelname)s %(message)s')
    LOG.info("BALSAMIC started with logger !!")


cli.add_command(install_command)
cli.add_command(run_command)
cli.add_command(config_command)
