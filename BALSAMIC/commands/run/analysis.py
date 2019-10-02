#!/usr/bin/env python
import os
import logging
import subprocess
import json
import click

# CLI commands and decorators
from BALSAMIC.utils.cli import createDir
from BALSAMIC.utils.cli import get_sbatchpy
from BALSAMIC.utils.cli import get_snakefile, SnakeMake
from BALSAMIC.utils.cli import get_config


LOG = logging.getLogger(__name__)


@click.command("analysis",
               short_help="Run the analysis on a provided sample config-file")
@click.option('-a',
              '--analysis-type',
              required=False,
              type=click.Choice(['qc', 'paired', 'single']),
              help='Type of analysis to run from input config file.\
              By default it will read from config file, but it will override config file \
              if it is set here.')
@click.option(
    '-S',
    '--snake-file',
    type=click.Path(),
    show_default=True,
    help='Input for a custom snakefile. WARNING: This is for internal testing,\
              and should not be used. Providing a snakefile supersedes analysis_type option.'
)
@click.option('-s',
              '--sample-config',
              required=True,
              type=click.Path(),
              help='Sample json config file.')
@click.option(
    '--run-mode',
    show_default=True,
    default='slurm',
    type=click.Choice(["local", "slurm"]),
    help='Run mode to use. By default SLURM will be used to run the analysis.\
              But local runner also available for local computing')
@click.option('-c',
              '--cluster-config',
              show_default=True,
              default=get_config('cluster'),
              type=click.Path(),
              help='SLURM config json file.')
@click.option('-l',
              '--log-file',
              type=click.Path(),
              help='Log file output for BALSAMIC.\
              This is raw log output from snakemake.')
@click.option(
    '-r',
    '--run-analysis',
    show_default=True,
    default=False,
    is_flag=True,
    help='By default balsamic run_analysis will run in dry run mode. \
              Raise thise flag to make the actual analysis')
@click.option('--qos',
              type=click.Choice(['low', 'normal', 'high']),
              show_default=True,
              default="low",
              help='QOS for sbatch jobs. Passed to ' + get_sbatchpy())
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
@click.option('--slurm-account', help='SLURM account to run jobs')
@click.option('--slurm-mail-user', help='SLURM mail user to send out email.')
@click.option(
    '--slurm-mail-type',
    type=click.Choice(
        ['NONE', 'BEGIN', 'END', 'FAIL', 'REQUEUE', 'ALL', 'TIME_LIMIT']),
    help='SLURM mail type to send out email.\
              This will be applied to all jobs and override snakemake settings.'
)
@click.pass_context
def analysis(context, snake_file, sample_config, run_mode, cluster_config,
             run_analysis, log_file, force_all, snakemake_opt, slurm_mail_type,
             slurm_mail_user, slurm_account, analysis_type, qos):
    """
    Runs BALSAMIC workflow on the provided sample's config file
    """

    if run_mode == 'slurm' and not run_analysis:
        LOG.info('Changing run-mode to local on dry-run')
        run_mode = 'local'

    if run_mode == 'slurm' and not slurm_account:
        LOG.info('slurm account is required for slurm run mode')
        raise click.Abort()

    sample_config_path = os.path.abspath(sample_config)

    with open(sample_config, 'r') as sample_fh:
        sample_config = json.load(sample_fh)

    logpath = sample_config['analysis']['log']
    scriptpath = sample_config['analysis']['script']
    resultpath = sample_config['analysis']['result']
    benchmarkpath = sample_config['analysis']['benchmark']
    case_name = sample_config['analysis']['case_id']
    sequencing_type = sample_config['analysis']['sequencing_type']

    if run_analysis:
        # if not dry run, then create (new) log/script directory
        for dirpath, dirnames, files in os.walk(logpath):
            if files:
                logpath = createDir(logpath, [])
                scriptpath = createDir(scriptpath, [])
                benchmarkpath = createDir(benchmarkpath, [])

    # Create result directory
    os.makedirs(resultpath, exist_ok=True)

    if not os.path.exists(logpath):
        os.makedirs(logpath, exist_ok=True)
        os.makedirs(scriptpath, exist_ok=True)
        os.makedirs(benchmarkpath, exist_ok=True)

    if not analysis_type:
        analysis_type = sample_config['analysis']['analysis_type']

    # Construct snakemake command to run workflow
    balsamic_run = SnakeMake()
    balsamic_run.case_name = case_name
    balsamic_run.working_dir = sample_config['analysis']['analysis_dir'] +  \
        case_name + '/BALSAMIC_run/'
    balsamic_run.snakefile = snake_file if snake_file else get_snakefile(
        analysis_type, sequencing_type)
    balsamic_run.configfile = sample_config_path
    balsamic_run.run_mode = run_mode
    balsamic_run.cluster_config = cluster_config
    balsamic_run.scheduler = get_sbatchpy()
    balsamic_run.log_path = logpath
    balsamic_run.script_path = scriptpath
    balsamic_run.result_path = resultpath
    balsamic_run.qos = qos
    balsamic_run.account = slurm_account
    if slurm_mail_type:
        balsamic_run.mail_type = slurm_mail_type
    balsamic_run.mail_user = slurm_mail_user
    balsamic_run.forceall = force_all
    balsamic_run.run_analysis = run_analysis
    balsamic_run.sm_opt = snakemake_opt

    try:
        subprocess.run(balsamic_run.build_cmd(), shell=True)
    except:
        raise 
