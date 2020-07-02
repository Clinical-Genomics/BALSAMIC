#!/usr/bin/env python
import os
import logging
import subprocess
import json
import click

# CLI commands and decorators
from BALSAMIC.utils.cli import createDir
from BALSAMIC.utils.cli import get_schedulerpy 
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
    default='cluster',
    type=click.Choice(["local", "cluster"]),
    help='Run mode to use. By default SLURM will be used to run the analysis.\
              But local runner also available for local computing')
@click.option('-c',
              '--cluster-config',
              show_default=True,
              default=get_config('cluster'),
              type=click.Path(),
              help='cluster config json file. (eg- SLURM, QSUB)')
@click.option('-l',
              '--log-file',
              type=click.Path(),
              help='Log file output for BALSAMIC.\
              This is raw log output from snakemake.')
@click.option('-p',
              '--profile',
              default="slurm",
              type=click.Choice(["slurm", "qsub"]),
              help="cluster profile to submit jobs")
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
              help='QOS for sbatch jobs. Passed to ' + get_schedulerpy())
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
@click.option('--account', '--slurm-account', '--qsub-account',
              help='cluster account to run jobs, ie: slurm_account')
@click.option('--mail-user', help='cluster mail user to send out email. e.g.: slurm_mail_user')
@click.option('--mail-type',
              type=click.Choice(['NONE', 'BEGIN', 'END', 'FAIL', 'REQUEUE', 'ALL', 'TIME_LIMIT']),
              help='cluster mail type to send out email. \
              This will be applied to all jobs and override snakemake settings.')
@click.pass_context
def analysis(context, snake_file, sample_config, run_mode, cluster_config,
             run_analysis, log_file, force_all, snakemake_opt, mail_type,
             mail_user, account, analysis_type, qos, profile):
    """
    Runs BALSAMIC workflow on the provided sample's config file
    """
    LOG.info(f"BALSAMIC started with log level {context.obj['loglevel']}.")

    if run_mode == 'cluster' and not run_analysis:
        LOG.info('Changing run-mode to local on dry-run')
        run_mode = 'local'

    if run_mode == 'cluster' and not account:
        LOG.info('slurm-account, qsub-account, or account is required for slurm run mode')
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
                sample_config['analysis']['benchmark'] = createDir(benchmarkpath, [])

    # Create result directory
    os.makedirs(resultpath, exist_ok=True)

    if not os.path.exists(logpath):
        os.makedirs(logpath, exist_ok=True)
        os.makedirs(scriptpath, exist_ok=True)
        os.makedirs(benchmarkpath, exist_ok=True)

    if not analysis_type:
        analysis_type = sample_config['analysis']['analysis_type']

    # Singularity bind path
    bind_path = list()
    bind_path.append(os.path.commonpath(sample_config['reference'].values()))
    if 'panel' in sample_config.keys():
        bind_path.append(sample_config['panel']['capture_kit'])
    bind_path.append(sample_config['analysis']['analysis_dir'])

    # Construct snakemake command to run workflow
    balsamic_run = SnakeMake()
    balsamic_run.case_name = case_name
    balsamic_run.working_dir = sample_config['analysis']['analysis_dir'] + "/" + \
        case_name + '/BALSAMIC_run/'
    balsamic_run.snakefile = snake_file if snake_file else get_snakefile(
        analysis_type, sequencing_type)
    balsamic_run.configfile = sample_config_path
    balsamic_run.run_mode = run_mode
    balsamic_run.cluster_config = cluster_config
    balsamic_run.scheduler = get_schedulerpy()
    balsamic_run.profile = profile
    balsamic_run.log_path = logpath
    balsamic_run.script_path = scriptpath
    balsamic_run.result_path = resultpath
    balsamic_run.qos = qos
    balsamic_run.account = account
    if mail_type:
        balsamic_run.mail_type = mail_type
    balsamic_run.mail_user = mail_user
    balsamic_run.forceall = force_all
    balsamic_run.run_analysis = run_analysis
    # Always use singularity
    balsamic_run.use_singularity = True
    balsamic_run.singularity_bind = bind_path
    balsamic_run.sm_opt = snakemake_opt

    try:
        subprocess.run(balsamic_run.build_cmd(), shell=True)#, check=True)
    except Exception as e:
        print(e)
        raise click.Abort()
