"""Balsamic commands init options."""
import click

OPTION_OUT_DIR = click.option(
    "-o",
    "--out-dir",
    required=True,
    type=click.Path(exists=True),
    help="Output directory for singularity containers and reference files",
)

OPTION_COSMIC_KEY = click.option(
    "-c",
    "--cosmic-key",
    required=False,
    type=click.STRING,
    help="Cosmic DB authentication key",
)

OPTION_CLUSTER_CONFIG = click.option(
    "--cluster-config",
    type=click.Path(),
    help="Cluster configuration JSON file path",
)
