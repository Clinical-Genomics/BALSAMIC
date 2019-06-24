#!/usr/bin/env python
import click

#from BALSAMIC.commands.config.sample import sample as sample_command
from BALSAMIC.commands.initiate.reference import reference as reference_command


@click.group()
@click.pass_context
def initiate(context):
    "create config files required for running the pipeline and reporting it"
    pass


#config.add_command(sample_command)
initiate.add_command(reference_command)
