#! /usr/bin/env python

import json
import subprocess
import logging
import click

from BALSAMIC.utils.cli import get_schedulerpy
from BALSAMIC.utils.cli import get_snakefile, SnakeMake, get_config

LOG = logging.getLogger(__name__)


@click.command('reference', short_help="Run the GenerateRef.smk workflow")
@click.option('-s',
              "--snakefile",
              default=get_snakefile('generate_ref'),
              help="snakefile for reference generation")
@click.option('-c',
              '--configfile',
              required=True,
              help="Config file to run the workflow")
@click.option('--run-mode',
              default='local',
              type=click.Choice(["local"]),
              help="Run mode to use. Only local supported for this.")
@click.option('--cluster-config',
              show_default=True,
              default=get_config('cluster'),
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
    help='By default balsamic run_analysis will run in dry run mode. \
              Raise thise flag to make the actual analysis')
@click.option(
    '-f',
    '--force-all',
    show_default=True,
    default=False,
    is_flag=True,
    help='Force run all analysis. This is same as snakemake --forceall')
@click.option('--snakemake-opt',
              multiple=True,
              help='Pass these options directly to snakemake')
@click.pass_context
def reference(context, snakefile, configfile, run_mode, cluster_config,
              log_file, run_analysis, force_all, snakemake_opt):
    """ Run generate reference workflow """
    LOG.info(f"BALSAMIC started with log level {context.obj['loglevel']}.")
    LOG.info("Reference generation workflow started")

    with open(configfile, "r") as config_fh:
        config = json.load(config_fh)

    # Singularity bind path
    bind_path = list()
    bind_path.append(config['output'])
    bind_path.append(config['conda_env_yaml'])
    bind_path.append(config['rule_directory'])

    # Construct snakemake command to run workflow
    balsamic_run = SnakeMake()
    balsamic_run.working_dir = config['output']
    balsamic_run.snakefile = snakefile
    balsamic_run.configfile = configfile
    balsamic_run.run_mode = run_mode
    balsamic_run.forceall = force_all
    balsamic_run.run_analysis = run_analysis
    balsamic_run.sm_opt = snakemake_opt
    # Always use singularity
    balsamic_run.use_singularity = True
    balsamic_run.singularity_bind = bind_path
    balsamic_run.sm_opt = snakemake_opt

    subprocess.run(balsamic_run.build_cmd(), shell=True)
