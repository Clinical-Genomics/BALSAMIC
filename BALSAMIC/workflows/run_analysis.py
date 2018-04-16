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
    '-d',
    '--dry-run',
    show_default=True,
    default=False,
    is_flag=True,
    help='Dryrun to check the output.')
@click.option(
    '-f',
    '--force-all',
    show_default=True,
    default=False,
    is_flag=True,
    help='Force run all analysis. This is same as snakemake --forceall')
@click.pass_context
def run_analysis(context, snake_file, config_sample, config_slurm, dry_run,
                 log_file, force_all):
    """

    Runs BALSAMIC workflow on the provided sample's config file
    
    """

    click.echo("Running the analysis")
    shellcmd = ["snakemake --immediate-submit -j 99 --notemp -np"]
    shellcmd.append("--jobname " + get_sample_name(config_sample) + ".{rulename}.{jobid}.py")
    shellcmd.append("--snakefile " + snake_file)
    shellcmd.append("--configfile " + config_sample)
    shellcmd.append("--cluster-config " + config_slurm)
    shellcmd.append("--cluster 'python3 runSbatch.py --sample_config " + config_sample + " {dependencies} '")
    print(" ".join(shellcmd))
    print(shellcmd)
    
    #subprocess.check_output(" ".join(shellcmd), shell=True, stderr=subprocess.STDOUT)
    subprocess.Popen(" ".join(shellcmd), shell=True)
