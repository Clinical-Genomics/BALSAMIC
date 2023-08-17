"""Balsamic commands init options."""
import click

from BALSAMIC.constants.cluster import ClusterConfigType
from BALSAMIC.utils.cli import get_config_path

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

OPTION_SNAKEFILE = click.option(
    "-s",
    "--snakefile",
    default=None,
    type=click.Path(),
    show_default=True,
    help="Snakefile for reference generation",
)

OPTION_CLUSTER_CONFIG = click.option(
    "--cluster-config",
    type=click.Path(),
    default=get_config_path(ClusterConfigType.CACHE),
    show_default=True,
    help="Cluster configuration JSON file path",
)
