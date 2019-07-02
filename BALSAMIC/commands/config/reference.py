#! /usr/bin/env python

import os
import click

from BALSAMIC.utils.cli import write_json


@click.command("reference", short_help="config workflow for generate reference")
@click.option("-o", "--outdir", required=True, help="output directory for ref files eg: reference/")
@click.option("-c", "--cosmic-key", help="cosmic db authentication key")
def reference(outdir, cosmic_key):
    """ Configure workflow for reference generation """
    config = dict()
    outdir = os.path.abspath(outdir)
    config_json = os.path.join(outdir, "config.json")

    config["output"] = outdir
    if cosmic_key:
        config["cosmic"] = cosmic_key

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    write_json(config, config_json)
