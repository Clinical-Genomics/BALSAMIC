#!/usr/bin/env python
import click

from BALSAMIC.commands.plugins.scout import scout as scout_command
from BALSAMIC.commands.plugins.cov_plot import target_cov_plot as target_cov_plot_command


@click.group()
@click.pass_context
def plugins(context):
    ''' Additional and helper utilities for third party applications '''
    pass


plugins.add_command(scout_command)
plugins.add_command(target_cov_plot_command)
