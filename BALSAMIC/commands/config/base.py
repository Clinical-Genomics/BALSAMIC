#!/usr/bin/env python
import click

from BALSAMIC.commands.config.case import case_config as case_command
from BALSAMIC.commands.config.pon import pon_config as pon_command

@click.group()
@click.pass_context
def config(context):
    "create config files required for running the pipeline."
    pass


config.add_command(case_command)
config.add_command(pon_command)
