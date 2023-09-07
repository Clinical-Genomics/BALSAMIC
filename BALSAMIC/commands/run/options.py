"""Balsamic run command options."""
import click


OPTION_DRAGEN = click.option(
    "--dragen",
    is_flag=True,
    default=False,
    help="Enable dragen variant caller",
)

OPTION_BENCHMARK = click.option(
    "--benchmark",
    default=False,
    is_flag=True,
    help="Profile slurm jobs. Make sure you have slurm profiler enabled in your HPC.",
)
