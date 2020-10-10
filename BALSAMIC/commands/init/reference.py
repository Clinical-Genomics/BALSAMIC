import os
import logging
import click
import graphviz
import snakemake
import subprocess
from pathlib import Path

from BALSAMIC.utils.cli import write_json, merge_json, CaptureStdout, get_snakefile, SnakeMake
from BALSAMIC import __version__ as balsamic_version

LOG = logging.getLogger(__name__)


@click.command("reference",
               short_help="config workflow for generate reference")
@click.option("-o",
              "--outdir",
              "--out-dir",
              required=True,
              help=("Output directory for ref files."
                    "This path will be used as base path for files"))
@click.option("-c",
              "--cosmic-key",
              required=True,
              help="cosmic db authentication key")
@click.option("-s",
              "--snakefile",
              default=get_snakefile('generate_ref'),
              type=click.Path(),
              show_default=True,
              help="snakefile for reference generation")
@click.option("-d",
              "--dagfile",
              default="generate_ref_worflow_graph",
              show_default=True,
              help="DAG file for overview")
@click.option("--singularity",
              type=click.Path(),
              required=True,
              help='Download singularity image for BALSAMIC')
@click.option("-g",
              "--genome-version",
              default="hg19",
              type=click.Choice(["hg19", "hg38"]),
              help=("Genome version to prepare reference. Path to genome"
                    "will be <outdir>/genome_version"))
@click.option(
    '-r',
    '--run-analysis',
    show_default=True,
    default=False,
    is_flag=True,
    help=("By default balsamic run_analysis will run in dry run mode."
          "Raise thise flag to make the actual analysis"))
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
@click.option('-q',
              '--quiet',
              default=False,
              is_flag=True,
              help=('Instruct snakemake to be quiet!'
                    'No output will be printed'))
@click.pass_context
def reference(context, outdir, cosmic_key, snakefile, dagfile, singularity,
              genome_version, run_analysis, force_all, quiet, snakemake_opt):
    """ Configure workflow for reference generation """

    LOG.info(f"BALSAMIC started with log level {context.obj['loglevel']}.")
    config_path = Path(__file__).parents[2] / "config"
    config_path = config_path.absolute()

    balsamic_env = config_path / "balsamic_env.yaml"
    rule_directory = Path(__file__).parents[2]

    install_config = dict()

    install_config["conda_env_yaml"] = balsamic_env.as_posix()
    install_config["rule_directory"] = rule_directory.as_posix() + "/"

    install_config["singularity"] = dict()
    install_config["singularity"]["image"] = Path(
        singularity).absolute().as_posix()

    config = dict()
    outdir = os.path.join(os.path.abspath(outdir), balsamic_version,
                          genome_version)
    config_json = os.path.join(outdir, "config.json")
    dagfile_path = os.path.join(outdir, dagfile)

    config["output"] = outdir
    if cosmic_key:
        config["cosmic_key"] = cosmic_key

    config["genome_version"] = genome_version

    config = merge_json(config, install_config)

    os.makedirs(outdir, exist_ok=True)

    write_json(config, config_json)
    LOG.info(
        f'Reference generation workflow configured successfully - {config_json}'
    )

    with CaptureStdout() as graph_dot:
        snakemake.snakemake(snakefile=snakefile,
                            dryrun=True,
                            configfiles=[config_json],
                            printrulegraph=True)

    graph_title = "_".join(
        ['BALSAMIC', balsamic_version, 'Generate reference'])
    graph_dot = "".join(graph_dot).replace(
        'snakemake_dag {',
        'BALSAMIC { label="' + graph_title + '";labelloc="t";')
    graph_obj = graphviz.Source(graph_dot,
                                filename=dagfile_path,
                                format="pdf",
                                engine="dot")

    try:
        graph_pdf = graph_obj.render()
        LOG.info(
            f'Reference workflow graph generated successfully - {graph_pdf}')
    except Exception:
        LOG.error(f'Reference workflow graph generation failed')
        raise click.Abort()

    LOG.info(f"BALSAMIC started with log level {context.obj['loglevel']}.")
    LOG.info("Reference generation workflow started")

    # Singularity bind path
    bind_path = list()
    bind_path.append(config['output'])
    bind_path.append(config['conda_env_yaml'])
    bind_path.append(config['rule_directory'])

    # Construct snakemake command to run workflow
    balsamic_run = SnakeMake()
    balsamic_run.working_dir = config['output']
    balsamic_run.snakefile = snakefile
    balsamic_run.configfile = config_json
    balsamic_run.run_mode = "local"
    balsamic_run.forceall = force_all
    balsamic_run.run_analysis = run_analysis
    balsamic_run.quiet = quiet
    balsamic_run.sm_opt = list(snakemake_opt).extend(["--cores", "1"])

    # Always use singularity
    balsamic_run.use_singularity = True
    balsamic_run.singularity_bind = bind_path

    subprocess.run(balsamic_run.build_cmd(), shell=True)
