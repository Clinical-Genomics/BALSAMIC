import os
import sys
import re
import logging
import subprocess
from pathlib import Path

import click
import graphviz
import snakemake

from BALSAMIC.utils.constants import (
    BIOINFO_TOOL_ENV,
    BALSAMIC_DOCKER_PATH,
    VALID_CONTAINER_CONDA_NAME,
)
from BALSAMIC.utils.cli import (
    write_json,
    merge_json,
    CaptureStdout,
    get_snakefile,
    SnakeMake,
    get_config,
    get_schedulerpy,
    job_id_dump_to_yaml,
)
from BALSAMIC import __version__ as balsamic_version

LOG = logging.getLogger(__name__)


@click.command(
    "init",
    short_help="Download matching version for container and build reference")
@click.option(
    "-o",
    "--outdir",
    "--out-dir",
    required=True,
    help=("Output directory for ref files."
          "This path will be used as base path for files"),
)
@click.option(
    "-v",
    "--container-version",
    show_default=True,
    default=balsamic_version,
    type=click.Choice(["develop", "master", balsamic_version]),
    help="Container for BALSAMIC version to download",
)
@click.option(
    "-f",
    "--force",
    show_default=True,
    default=False,
    is_flag=True,
    help="Force re-downloading all containers",
)
@click.option("-c",
              "--cosmic-key",
              required=True,
              help="cosmic db authentication key")
@click.option(
    "-s",
    "--snakefile",
    default=get_snakefile("generate_ref"),
    type=click.Path(),
    show_default=True,
    help="snakefile for reference generation",
)
@click.option(
    "-d",
    "--dagfile",
    default="generate_ref_worflow_graph",
    show_default=True,
    help="DAG file for overview",
)
@click.option(
    "-g",
    "--genome-version",
    default="hg19",
    type=click.Choice(["hg19", "hg38"]),
    help=("Genome version to prepare reference. Path to genome"
          "will be <outdir>/genome_version"),
)
@click.option(
    "-r",
    "--run-analysis",
    show_default=True,
    default=False,
    is_flag=True,
    help=("By default balsamic run_analysis will run in dry run mode."
          "Raise this flag to make the actual analysis"),
)
@click.option(
    "--run-mode",
    show_default=True,
    default="cluster",
    type=click.Choice(["local", "cluster"]),
    help=
    ("Run mode to use. By default SLURM will be used to generate the balsamic_cache"
     "Alternatively, option for local computing"),
)
@click.option(
    "-c",
    "--cluster-config",
    show_default=True,
    default=get_config("cluster"),
    type=click.Path(),
    help="cluster config json file. (eg- SLURM, QSUB)",
)
@click.option(
    "-p",
    "--profile",
    default="slurm",
    type=click.Choice(["slurm", "qsub"]),
    help="cluster profile to submit jobs",
)
@click.option(
    "--account",
    "--slurm-account",
    "--qsub-account",
    help="cluster account to run jobs, ie: slurm_account",
)
@click.option(
    "--qos",
    type=click.Choice(["low", "normal", "high"]),
    show_default=True,
    default="low",
    help="QOS for sbatch jobs. Passed to " + get_schedulerpy(),
)
@click.option("--mail-user",
              help="cluster mail user to send out email. e.g.: slurm_mail_user"
              )
@click.option(
    "--mail-type",
    type=click.Choice([
        "NONE",
        "BEGIN",
        "END",
        "FAIL",
        "REQUEUE",
        "ALL",
        "TIME_LIMIT",
    ]),
    help=("cluster mail type to send out email. This will "
          "be applied to all jobs and override snakemake settings."),
)
@click.option(
    "-f",
    "--force-all",
    show_default=True,
    default=False,
    is_flag=True,
    help="Force run all analysis. This is same as snakemake --forceall",
)
@click.option("--snakemake-opt",
              multiple=True,
              help="Pass these options directly to snakemake")
@click.option(
    "-q",
    "--quiet",
    default=False,
    is_flag=True,
    help="Instruct snakemake to be quiet! No output will be printed",
)
@click.pass_context
def initialize(
    context,
    outdir,
    container_version,
    force,
    cosmic_key,
    snakefile,
    dagfile,
    genome_version,
    run_analysis,
    run_mode,
    cluster_config,
    account,
    qos,
    profile,
    mail_user,
    mail_type,
    force_all,
    quiet,
    snakemake_opt,
):
    """
    Initialize various resources after first installation.
    - Pull container(s) for BALSAMIC according to matching version
    - Download and build a reference
    """
    LOG.info("BALSAMIC started with log level %s" % context.obj["loglevel"])

    if run_mode == "cluster" and not run_analysis:
        LOG.info("Changing run-mode to local on dry-run")
        run_mode = "local"

    if run_mode == "cluster" and not account:
        LOG.info(
            "slurm-account, qsub-account, or account is required for slurm run mode"
        )
        raise click.Abort()


#    sbatch_template = ("#!/bin/bash -l\n"
#                       "#SBATCH --account={account_type}\n"
#                       "#SBATCH -N 1 -c 36 \n"
#                       "#SBATCH -t 5:00:00\n"
#                       "#SBATCH --qos={qos}\n"
#                       "#SBATCH -J balsamic_cache\n")

# resolve outdir to absolute path
    outdir = Path(outdir).resolve()
    #scripts_dir = Path(outdir, "scripts")
    #Path(scripts_dir).mkdir(parents=True, exist_ok=True)

    container_outdir = Path(outdir, balsamic_version, "containers")
    pattern = re.compile(r"^(\d+\.)?(\d+\.)?(\*|\d+)$")
    if pattern.findall(container_version):
        docker_image_base_name = "release_v{}".format(container_version)
    else:
        docker_image_base_name = container_version

    for image_suffix in VALID_CONTAINER_CONDA_NAME:
        container_stub_url = "{}:{}-{}".format(BALSAMIC_DOCKER_PATH,
                                               docker_image_base_name,
                                               image_suffix)
        # Pull container
        LOG.info("Singularity image source: {}".format(container_stub_url))

        # Set container name according to above docker image name
        Path(container_outdir).mkdir(parents=True, exist_ok=True)
        image_name = Path(container_outdir,
                          "{}.sif".format(image_suffix)).as_posix()
        LOG.info("Image will be downloaded to {}".format(image_name))
        LOG.info("Starting download. This process can take some time...")

        cmd = ["singularity", "pull", "--name", f"{image_name}"]
        if force:
            cmd.append("--force")
        cmd.append(container_stub_url)

        LOG.info("The following command will run: {}".format(" ".join(cmd)))
        if run_analysis:
            subprocess.run(" ".join(cmd), shell=True)

    # job_fl = os.path.join(
    #    scripts_dir, "balsamic_container_{}.sh".format(image_suffix)
    # )
    # template = sbatch_template + "\n" + " ".join(cmd)
    # with open(job_fl, "w") as jb_fl:
    #    jb_fl.write(template.format(account_type=account, qos=qos))
    # LOG.info("The following command will run: {}".format(" ".join(cmd)))
    # if run_analysis and run_mode == "cluster":
    #    subprocess.check_call(["sbatch", job_fl])
    # elif run_analysis and run_mode == "local":
    #    subprocess.run(" ".join(cmd), shell=True)

    config_path = Path(__file__).parents[2] / "config"
    config_path = config_path.absolute()

    rule_directory = Path(__file__).parents[2]

    config_dict = dict()
    config_dict["bioinfo_tools"] = BIOINFO_TOOL_ENV
    config_dict["rule_directory"] = rule_directory.as_posix() + "/"
    config_dict["singularity"] = dict()
    config_dict["singularity"]["image"] = container_outdir.as_posix()

    reference_outdir = Path(outdir, balsamic_version, genome_version)
    Path(reference_outdir).mkdir(parents=True, exist_ok=True)
    config_json = Path(reference_outdir, "config.json").as_posix()
    dagfile_path = Path(reference_outdir, dagfile).as_posix()
    logpath = Path(reference_outdir, "logs")
    Path(logpath).mkdir(parents=True, exist_ok=True)
    scriptpath = Path(reference_outdir, "scripts")
    Path(scriptpath).mkdir(parents=True, exist_ok=True)

    config_dict["output"] = reference_outdir.as_posix()
    if cosmic_key:
        config_dict["cosmic_key"] = cosmic_key

    config_dict["genome_version"] = genome_version

    write_json(config_dict, config_json)
    LOG.info("Reference generation workflow configured successfully - %s" %
             config_json)

    with CaptureStdout() as graph_dot:
        snakemake.snakemake(
            snakefile=snakefile,
            dryrun=True,
            configfiles=[config_json],
            printrulegraph=True,
        )

    graph_title = "_".join(
        ["BALSAMIC", balsamic_version, "Generate reference"])
    graph_dot = "".join(graph_dot).replace(
        "snakemake_dag {",
        'BALSAMIC { label="' + graph_title + '";labelloc="t";')
    graph_obj = graphviz.Source(graph_dot,
                                filename=dagfile_path,
                                format="pdf",
                                engine="dot")

    try:
        graph_pdf = graph_obj.render()
        LOG.info("Reference workflow graph generated successfully - %s " %
                 graph_pdf)
    except Exception:
        LOG.error("Reference workflow graph generation failed")
        raise click.Abort()

    LOG.info("Reference generation workflow started")

    # Singularity bind path
    bind_path = list()
    bind_path.append(config_dict["output"])
    bind_path.append(config_dict["rule_directory"])

    # Construct snakemake command to run workflow
    balsamic_run = SnakeMake()
    balsamic_run.working_dir = config_dict["output"]
    balsamic_run.snakefile = snakefile
    balsamic_run.configfile = config_json
    balsamic_run.run_mode = run_mode
    balsamic_run.forceall = force_all
    balsamic_run.run_analysis = run_analysis
    balsamic_run.cluster_config = cluster_config
    balsamic_run.scheduler = get_schedulerpy()
    balsamic_run.profile = profile
    balsamic_run.account = account
    balsamic_run.qos = qos
    balsamic_run.log_path = logpath
    balsamic_run.script_path = scriptpath
    balsamic_run.quiet = quiet
    balsamic_run.sm_opt = snakemake_opt  # + ["--cores", "1"]

    # Always use singularity
    balsamic_run.use_singularity = True
    balsamic_run.singularity_bind = bind_path

    cmd = sys.executable + " -m " + balsamic_run.build_cmd()
    #subprocess.run(balsamic_run.build_cmd(), shell=True)
    #subprocess.run(cmd, shell=True)

    if run_analysis and run_mode == "cluster":
        jobid_dump = os.path.join(logpath, "reference_generation.sacct")
        jobid_yaml = os.path.join(reference_outdir, profile + "_jobids.yaml")
        job_id_dump_to_yaml(jobid_dump, jobid_yaml, "reference_generation")
    elif run_analysis and run_mode == "local":
        subprocess.run(cmd, shell=True)
