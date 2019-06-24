#!/usr/bin/env python
import click

from BALSAMIC import __version__ as bv


@click.command("reference",
               short_help="Downloads and prepares reference files for BALSAMIC"
               )
@click.option(
    "-o",
    "--output-config",
    required=False,
    help="A file containing the output config reference files."
    "This file ,as reference.json, is required for running BALSAMIC.",
)
@click.option(
    "-d",
    "--output-dir",
    required=True,
    type=click.Path(),
    help="Output directory to store reference files.",
)
@click.pass_context
def reference(
        context,
        output_config,
        output_dir,
):
    """
    Downloads and creates reference files for BALSAMIC. Requires a compelete BALSAMIC installation.
    """
    pass
