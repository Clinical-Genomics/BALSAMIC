"""Balsamic CLI."""
import logging

import click
import coloredlogs

from BALSAMIC import __version__ as balsamic_version
from BALSAMIC.commands.config.base import config as config_command
from BALSAMIC.commands.init.base import initialize as init_command
from BALSAMIC.commands.plugins.base import plugins as plugins_command
from BALSAMIC.commands.report.base import report as report_command
from BALSAMIC.commands.run.base import run as run_command
from BALSAMIC.constants.constants import LogLevel, LOG_LEVELS
from BALSAMIC.utils.cli import add_doc as doc

LOG = logging.getLogger(__name__)


@click.group()
@click.option(
    "--log-level",
    default=LogLevel.INFO.value,
    type=click.Choice(LOG_LEVELS),
    help="Logging level in terms of urgency",
    show_default=True,
)
@click.version_option(version=balsamic_version)
@click.pass_context
@doc(
    f"Balsamic {balsamic_version}: Bioinformatic Analysis Pipeline for Somatic Mutations in Cancer"
)
def cli(context, log_level):
    coloredlogs.DEFAULT_FIELD_STYLES = {
        "asctime": {"color": "green"},
        "hostname": {"color": "magenta"},
        "levelname": {"color": "yellow", "bold": True},
        "programname": {"color": "cyan"},
        "name": {"color": "blue"},
    }
    coloredlogs.install(
        level=log_level,
        fmt="%(programname)s %(hostname)s %(asctime)s %(name)s pid:%(process)d [%(levelname)s] %(message)s",
    )
    LOG.info(f"Running BALSAMIC version {balsamic_version}")
    context.obj["log_level"] = log_level


cli.add_command(run_command)
cli.add_command(report_command)
cli.add_command(config_command)
cli.add_command(plugins_command)
cli.add_command(init_command)
