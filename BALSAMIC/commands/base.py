"""Balsamic CLI."""
import logging

import click
import coloredlogs

from BALSAMIC import __version__ as balsamic_version
from BALSAMIC.commands.config.base import config as config_command
from BALSAMIC.commands.init.base import initialize as init_command
from BALSAMIC.commands.options import OPTION_LOG_LEVEL
from BALSAMIC.commands.report.base import report as report_command
from BALSAMIC.commands.run.base import run as run_command
from BALSAMIC.constants.constants import LogLevel
from BALSAMIC.utils.cli import add_doc as doc

LOG = logging.getLogger(__name__)


@click.group()
@OPTION_LOG_LEVEL
@click.version_option(version=balsamic_version)
@click.pass_context
@doc(
    f"Balsamic {balsamic_version}: Bioinformatic Analysis Pipeline for Somatic Mutations in Cancer"
)
def cli(context: click.Context, log_level: LogLevel):
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
    context.obj = {"log_level": log_level}


cli.add_command(run_command)
cli.add_command(report_command)
cli.add_command(config_command)
cli.add_command(init_command)
