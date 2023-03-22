import click

from BALSAMIC.constants.analysis import ConfigType
from BALSAMIC.constants.cache import ContainerVersion
from BALSAMIC.utils.cli import get_config_path

OPTION_OUT_DIR = click.option(
    "-o",
    "--out-dir",
    "--outdir",
    required=True,
    type=click.Path(exists=True),
    help="Output directory for singularity containers and reference files",
)

OPTION_CONTAINER_VERSION = click.option(
    "-v",
    "--container-version",
    show_default=True,
    default=ContainerVersion.RELEASE,
    type=click.Choice([ContainerVersion.DEVELOP, ContainerVersion.RELEASE]),
    help="Container version to be downloaded",
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
    default=get_config_path(ConfigType.CLUSTER_REFERENCE),
    show_default=True,
    help="Cluster configuration JSON file path",
)
