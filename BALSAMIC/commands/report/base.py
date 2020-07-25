#!/usr/bin/env python
import click

from BALSAMIC.commands.report.deliver import deliver as deliver_command
from BALSAMIC.commands.report.status import status as status_command


@click.group()
@click.pass_context
def report(context):
    ''' Various command to create report, check status, and prepare delivery files '''
    pass


report.add_command(deliver_command)
report.add_command(status_command)
