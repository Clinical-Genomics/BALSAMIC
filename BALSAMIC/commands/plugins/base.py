#!/usr/bin/env python
import click

from BALSAMIC.commands.plugins.housekeeper import housekeeper as housekeeper_command
from BALSAMIC.commands.plugins.scout import scout as scout_command


@click.group()
@click.pass_context
def plugins(context):
    ''' Additional and helper utilities for third party applications '''
    pass


plugins.add_command(housekeeper_command)
plugins.add_command(scout_command)
