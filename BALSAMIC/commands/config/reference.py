#! /usr/bin/env python

import os
import subprocess
import click

from BALSAMIC.utils.cli import write_json
from BALSAMIC.utils.cli import get_snakefile, get_config
from BALSAMIC.utils.cli import SnakeMake
from BALSAMIC.commands.config.sample import merge_json
from BALSAMIC import __version__ as bv


@click.command("reference",
               short_help="config workflow for generate reference")
@click.option("-o",
              "--outdir",
              required=True,
              help="output directory for ref files eg: reference")
@click.option("-i",
              "--install-config",
              type=click.Path(),
              help="install config file.")
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
              default="generate_ref_dag.pdf",
              show_default=True,
              help="DAG file for overview")
@click.option("--singularity",
              is_flag=True,
              help="To run the workflow on container")
def reference(outdir, cosmic_key, snakefile, dagfile, singularity,
              install_config):
    """ Configure workflow for reference generation """
    if not install_config:
        try:
            install_config = get_config("install")
        except:
            context.Abort() 

    config = dict()
    outdir = os.path.abspath(outdir)
    config_json = os.path.join(outdir, "config.json")
    dagfile_path = os.path.join(outdir, dagfile)

    config["output"] = outdir
    if cosmic_key:
        config["cosmic_key"] = cosmic_key

    config = merge_json(config, install_config)

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    write_json(config, config_json)

    # configure snakemake cmd
    config_reference = SnakeMake()
    config_reference.working_dir = outdir
    config_reference.snakefile = snakefile
    config_reference.configfile = config_json

    # To create DAG file
    shell_cmd = [
        '"--rulegraph"', "|", "sed", '"s/digraph', 'snakemake_dag',
        '{/digraph', 'BALSAMIC', '{', 'labelloc=\\"t\\"\;', 'label=\\"Title:',
        'Reference', 'generation', 'for', 'BALSAMIC', bv, 'workflow',
        '\\"\;/g"', '|', 'dot', '-Tpdf', '1>', dagfile_path
    ]

    cmd = config_reference.build_cmd() + " " + " ".join(shell_cmd)
    click.echo("Creating workflow DAG image: %s" % dagfile_path)
    subprocess.run(cmd, shell=True)
