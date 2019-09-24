#! /usr/bin/env python

import os
import logging
import click
import graphviz
import snakemake

from BALSAMIC.utils.cli import write_json
from BALSAMIC.utils.cli import get_snakefile, get_config
from BALSAMIC.utils.cli import CaptureStdout
from BALSAMIC.commands.config.case import merge_json
from BALSAMIC import __version__ as bv

LOG = logging.getLogger(__name__)


@click.command("reference",
               short_help="config workflow for generate reference")
@click.option("-o",
              "--outdir",
              required=True,
              help="output directory for ref files eg: reference")
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
def reference(outdir, cosmic_key, snakefile, dagfile):
    """ Configure workflow for reference generation """

    install_config = get_config("install")

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

    with CaptureStdout() as graph_dot:
        snakemake.snakemake(snakefile=snakefile,
                            dryrun=True,
                            configfile=config_json,
                            printrulegraph=True)

    graph_title = "_".join(['BALSAMIC', bv, 'Generate reference'])
    graph_dot = "".join(graph_dot).replace(
        'snakemake_dag {',
        'BALSAMIC { label="' + graph_title + '";labelloc="t";')
    graph_obj = graphviz.Source(graph_dot,
                                filename=dagfile_path,
                                format="pdf",
                                engine="dot")

    if graph_obj.render():
        LOG.info(f'Reference generation workflow configured successfully - {outdir}')
    else:
        LOG.error(f'Snakemake DAG graph generation failed - {dagfile_path}')
        click.Abort()
