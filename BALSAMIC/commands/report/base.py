"""Balsamic report CLI."""
import click

from BALSAMIC.commands.report.deliver import deliver as deliver_command

@click.group()
@click.pass_context
def report(context: click.Context):
    """Command to generate delivery files and check analysis status."""
    pass

report.add_command(deliver_command)

