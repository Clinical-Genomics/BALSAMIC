#!/usr/bin/env python
import click

from BALSAMIC.commands.plugins.delivery import delivery as delivery_command

@click.group()
@click.pass_context
def plugins(context):
    ''' Additional and helper utilities for third party applications '''
    pass


plugins.add_command(delivery_command)
