#!/usr/bin/env python
import os
import subprocess
import click
import logging
import sys
import json


def get_sample_name(json_in):

    json_in = json.load(open(json_in))

    return json_in["analysis"]["sample_id"]


@click.command(
    "run_analysis", short_help="Run BALSAMIC on a provided config file")
@click.option(
    '-S',
    '--snake-file',
    required=True,
    type=click.Path(),
    help='Snakefile required for snakemake to function.')
@click.option(
    '-c',
    '--config-sample',
    required=True,
    type=click.Path(),
    help='Sample json config file.')
@click.option(
    '-s',
    '--config-slurm',
    required=True,
    type=click.Path(),
    help='Slurm config json file.')
@click.option(
    '-l',
    '--log-file',
    type=click.Path(),
    help='Log file output for BALSAMIC. This is raw log output from snakemake.'
)
@click.option(
    '-r',
    '--run-analysis',
    show_default=True,
    default=False,
    is_flag=True,
    help='By default balsamic run_analysis will run in dry run mode. Raise thise flag to make the actual analysis')
@click.option(
    '-f',
    '--force-all',
    show_default=True,
    default=False,
    is_flag=True,
    help='Force run all analysis. This is same as snakemake --forceall')
@click.pass_context
def run_analysis(
        context,
        snake_file,
        config_sample,
        config_slurm,
        run_analysis,
        log_file,
        force_all):
    """

    Runs BALSAMIC workflow on the provided sample's config file

    """

    click.echo("Running the analysis")
    shellcmd = ["snakemake --immediate-submit -j 99 --notemp -p"]
    shellcmd.append(
        "--jobname " +
        get_sample_name(config_sample) +
        ".{rulename}.{jobid}.sh")
    shellcmd.append("--snakefile " + snake_file)
    shellcmd.append("--configfile " + config_sample)
    shellcmd.append("--cluster-config " + config_slurm)
    shellcmd.append(
        "--cluster 'python3 runSbatch.py --sample-config " +
        config_sample +
        " {dependencies} '")

    if not run_analysis:
        shellcmd.append("--dryrun")

    if force_all:
        shellcmd.append("--forceall")

    subprocess.run(" ".join(shellcmd), shell=True)
