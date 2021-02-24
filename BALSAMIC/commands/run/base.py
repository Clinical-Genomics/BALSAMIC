import click

from BALSAMIC.commands.run.analysis import analysis as run_analysis_cmd


@click.group()
@click.pass_context
def run(context):
    "Run BALSAMIC on a provided config file"
    pass


run.add_command(run_analysis_cmd)
