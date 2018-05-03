#!/usr/bin/env python
import os
import subprocess
import click
import logging
import sys
import json


def get_sample_name(json_in):
    """
    Get sample name from input json file
    """

    try:
        with open(os.path.abspath(json_in), "r") as fn:
            sample_json = json.load(fn)
            sample_name = sample_json["analysis"]["sample_id"]
    except OSError:
        print("Couldn't load json file or file path is not absolute")

    return sample_name


def get_sbatchpy():
    """
    Returns a string path for runSbatch.py
    """

    try:
        sbatch = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), 'runSbatch.py')
    except OSError:
        print("Couldn't locate sbatch submitter.")

    return sbatch


def get_snakefile():
    """
    Return a string path for variant calling snakefile.
    """

    try:
        snakefile = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), 'VariantCalling')
    except OSError:
        print("Couldn't locate variant calling snakefile.")

    return snakefile


def get_cluster_config():
    """
    Returns a string path for cluster config file.
    """

    try:
        cluster_config = os.path.abspath(
            os.path.join(
                os.path.dirname(os.path.abspath(__file__)), os.pardir,
                'config/cluster.json'))
    except OSError:
        print("Couldn't locate cluster config file.")

    return cluster_config


@click.command(
    "run_analysis", short_help="Run BALSAMIC on a provided config file")
@click.option(
    '-S',
    '--snake-file',
    type=click.Path(),
    show_default=True,
    default=get_snakefile(),
    help='Snakefile required for snakemake to function.')
@click.option(
    '-s',
    '--sample-config',
    required=True,
    type=click.Path(),
    help='Sample json config file.')
@click.option(
    '-c',
    '--cluster-config',
    show_default=True,
    default=get_cluster_config(),
    type=click.Path(),
    help='SLURM config json file.')
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
    help=
    'By default balsamic run_analysis will run in dry run mode. Raise thise flag to make the actual analysis'
)
@click.option(
    '-f',
    '--force-all',
    show_default=True,
    default=False,
    is_flag=True,
    help='Force run all analysis. This is same as snakemake --forceall')
@click.option(
    '--snakemake-opt', help='Pass these options directly to snakemake')
@click.pass_context
def run_analysis(context, snake_file, sample_config, cluster_config,
                 run_analysis, log_file, force_all, snakemake_opt):
    """

    Runs BALSAMIC workflow on the provided sample's config file

    """

    shellcmd = ["snakemake --immediate-submit -j 99 --notemp -p"]
    shellcmd.append("--jobname " + get_sample_name(sample_config) +
                    ".{rulename}.{jobid}.sh")
    shellcmd.append("--snakefile " + snake_file)
    shellcmd.append("--configfile " + sample_config)
    shellcmd.append("--cluster-config " + cluster_config)
    shellcmd.append("--cluster 'python3 " + get_sbatchpy() +
                    " --sample-config " + sample_config + " {dependencies} '")

    if not run_analysis:
        shellcmd.append("--dryrun")

    if force_all:
        shellcmd.append("--forceall")

    if snakemake_opt:
        shellcmd.append(" " + snakemake_opt)

    print(" ".join(shellcmd))


#    subprocess.run(" ".join(shellcmd), shell=True)
