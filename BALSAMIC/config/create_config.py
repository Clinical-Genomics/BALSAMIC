#!/usr/bin/env python
import os
import subprocess
import hashlib
import yaml
import click
import logging
import sys
import json


@click.command(
    "create_config",
    short_help="Create a sample config file from input sample data")
@click.option(
    '-s',
    '--sample-config',
    required=True,
    type=click.Path(),
    help='Input sample config file.')
@click.option(
    '-r',
    '--reference-config',
    required=True,
    type=click.Path(),
    help='Reference config file.')
@click.option(
    '-o',
    '--output-config',
    required=True,
    type=click.Path(),
    help='Output a json config file ready to be imported for run-analysis'
)
@click.pass_context
def create_config(context, sample_config, reference_config, output_config):
    """
    Prepares a config file for balsamic run_analysis.

    """

    click.echo("Reading sample config file %s" % sample_config)
    click.echo("Reading reference config file %s" % reference_config)
    click.echo("Writing output config file %s" % output_config)
