#!/usr/bin/env python
import os
import subprocess
import hashlib
import yaml
import click
import logging
import sys
import json

# Subcommands
from BALSAMIC.install import install_env as install_command
from BALSAMIC.workflows import run_analysis as analysis_command
from BALSAMIC.config import config as config_command
from BALSAMIC.dev_report import report as report_command

# CLI commands and decorators
from BALSAMIC.tools.cli_utils import add_doc as doc

# Get version
from BALSAMIC import __version__

@click.group()
@click.version_option(version=__version__)
@click.pass_context

@doc("BALSAMIC {version}: Bioinformatic Analysis pipeLine for SomAtic MutatIons in Cancer".format(version=__version__))
def cli(context):
    pass

cli.add_command(install_command)
cli.add_command(analysis_command)
cli.add_command(config_command)
cli.add_command(report_command)
