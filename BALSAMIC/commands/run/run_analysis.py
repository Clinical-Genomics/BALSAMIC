#!/usr/bin/env python
import os
import subprocess
import click
import logging
import sys
import json
from pathlib import Path

# CLI commands and decorators
from BALSAMIC.tools.cli_utils import createDir
from BALSAMIC.tools.cli_utils import get_sample_name
from BALSAMIC.tools.cli_utils import get_analysis_dir
from BALSAMIC.commands.config.sample import get_config


def get_analysis_type(json_in):
    """
    Get analysis type from input json file: single or paired
    """

    try:
        with open(os.path.abspath(json_in), "r") as fn:
            sample_json = json.load(fn)
            analysis_type = sample_json["analysis"]["analysis_type"]
    except OSError:
        print("Couldn't load json file or file path is not absolute")

    return analysis_type


def get_sbatchpy():
    """
    Returns a string path for runSbatch.py
    """

    try:
        p = Path(__file__).parents[0]
        sbatch = str(Path(p, 'runSbatch.py'))
    except OSError:
        print("Couldn't locate sbatch submitter.")

    return sbatch


def get_snakefile(analysis_type):
    """
    Return a string path for variant calling snakefile.
    """

    try:
        p = Path(__file__).parents[2]
        if analysis_type == "paired":
            snakefile = Path(p, 'workflows', 'VariantCalling_paired')
        elif analysis_type == "single":
            snakefile = Path(p, 'workflows', 'VariantCalling_single')
        elif analysis_type == "qc":
            snakefile = Path(p, 'workflows', 'Alignment')
        elif analysis_type == "paired_umi":
            snakefile = Path(p, 'workflows', 'VariantCalling_paired_umi')
        elif analysis_type == "single_umi":
            snakefile = Path(p, 'workflows', 'VariantCalling_single_umi')
        else:
            raise ValueError("analysis_type should be single or paired")

    except OSError:
        print("Couldn't locate variant calling snakefile.")

    return str(snakefile)


@click.command("run", short_help="Run BALSAMIC on a provided config file")
@click.option(
    '-a',
    '--analysis-type',
    required=False,
    type=click.Choice(['qc', 'paired', 'single', 'paired_umi']),
    help=
    'Type of analysis to run from input config file. By default it will read from config file, but it will override config file if it is set here.'
)
@click.option(
    '-S',
    '--snake-file',
    type=click.Path(),
    show_default=True,
    help=
    'Input for a custom snakefile. WARNING: This is for internal testing, and should not be used. Providing a snakefile supersedes analysis_type option.'
)
@click.option('-s',
              '--sample-config',
              required=True,
              type=click.Path(),
              help='Sample json config file.')
@click.option(
    '-R',
    '--runner',
    show_default=True,
    default='local',
    type=click.Choice(["local", "slurm"]),
    help=
    'Runner to use. By default local shell runner will be used to run the analysis. slurm runner also available for cluster computing'
)
@click.option('-c',
              '--cluster-config',
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
    help=
    'By default balsamic run_analysis will run in dry run mode. Raise thise flag to make the actual analysis'
)
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
              help='Pass these options directly to snakemake')
@click.pass_context
def run_analysis(context, snake_file, sample_config, runner, cluster_config,
                 run_analysis, log_file, force_all, snakemake_opt,
                 analysis_type, qos):
    """

    Runs BALSAMIC workflow on the provided sample's config file

    """

    # get result, log, script directories
    logpath = os.path.join(get_analysis_dir(sample_config, "analysis_dir"),
                           get_sample_name(sample_config),
                           get_analysis_dir(sample_config, "log"))
    scriptpath = os.path.join(get_analysis_dir(sample_config, "analysis_dir"),
                              get_sample_name(sample_config),
                              get_analysis_dir(sample_config, "script"))
    resultpath = os.path.join(get_analysis_dir(sample_config, "analysis_dir"),
                              get_sample_name(sample_config),
                              get_analysis_dir(sample_config, "result"))

    # Create result directory
    os.makedirs(resultpath, exist_ok=True)
    if not os.path.exists(logpath):
        os.makedirs(logpath, exist_ok=True)
        os.makedirs(scriptpath, exist_ok=True)

    shellcmd = ["snakemake --notemp -p"]

    shellcmd.append(
        "--directory " +
        os.path.join(get_analysis_dir(sample_config, "analysis_dir"),
                     get_sample_name(sample_config), "BALSAMIC_run"))

    if not analysis_type:
        analysis_type = get_analysis_type(sample_config)

    snakefile = snake_file if snake_file else get_snakefile(analysis_type)

    shellcmd.append("--snakefile " + snakefile)

    shellcmd.append("--configfile " + sample_config)

    if not run_analysis:
        shellcmd.append("--dryrun")

    if runner == 'slurm':
        shellcmd.append(" --immediate-submit -j 300 ")
        shellcmd.append("--jobname " + get_sample_name(sample_config) +
                        ".{rulename}.{jobid}.sh")
        shellcmd.append("--cluster-config " + cluster_config)
        shellcmd.append("--cluster 'python3 " + get_sbatchpy() +
                        " --sample-config " + os.path.abspath(sample_config) +
                        " --qos " + qos + " --dir-log " + logpath +
                        " --dir-script " + scriptpath + " --dir-result " +
                        resultpath + " {dependencies} '")

    if run_analysis:
        # if not dry run, then create (new) log/script directory
        for dirpath, dirnames, files in os.walk(logpath):
            if files:
                logpath = createDir(logpath, [])
                scriptpath = createDir(scriptpath, [])

    if force_all:
        shellcmd.append("--forceall")

    if snakemake_opt:
        shellcmd.append(" " + snakemake_opt)

    subprocess.run(" ".join(shellcmd), shell=True)
