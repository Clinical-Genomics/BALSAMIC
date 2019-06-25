#!/usr/bin/env python
import click

from BALSAMIC.commands.config.sample import sample as sample_command

@click.group()
@click.pass_context


def config(context):
    "create config files required for running the pipeline."
    pass


config.add_command(sample_command)
