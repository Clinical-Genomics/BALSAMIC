#!/usr/bin/env python
import os
import subprocess
import hashlib
import yaml
import click
import logging
import sys
import json

from BALSAMIC.install import install_conda as install_command
from BALSAMIC.workflows import run_analysis as run_analysis


@click.group()
@click.option(
    '-v',
    '--verbose',
    show_default=True,
    is_flag=True,
    help='Verbose output *Not implented yet')
@click.pass_context
def cli(context, verbose):
    """
    BALSAMIC: Bioinformatic Analysis pipeLine for SomAtic MutatIons in Cancer

    """
    pass


cli.add_command(install_command)
cli.add_command(run_analysis)
