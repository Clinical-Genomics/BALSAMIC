import click

from BALSAMIC.commands.init.reference import reference as reference_command
from BALSAMIC.commands.init.container import container as container_command


@click.group("init")
@click.pass_context
@click.option("-o",
              "--outdir",
              "--out-dir",
              required=True,
              help=("Output directory for ref files."
                    "This path will be used as base path for files"))
def initialize(context, outdir):
    "Initialize various resources after first installation."
    pass


initialize.add_command(reference_command)
initialize.add_command(container_command)
