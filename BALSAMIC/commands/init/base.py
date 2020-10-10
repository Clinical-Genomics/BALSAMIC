import click

from BALSAMIC.commands.init.reference import reference as reference_command


@click.group("init")
@click.pass_context
def initialize(context):
    "Initialize various resources after first installation."
    pass


initialize.add_command(reference_command)
