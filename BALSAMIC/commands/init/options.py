import click
from BALSAMIC import __version__ as balsamic_version

from BALSAMIC.constants.cache import ContainerVersion

OPTION_OUT_DIR = click.option(
    "-o",
    "--out-dir",
    "--outdir",
    required=True,
    type=click.Path(exists=True),
    help="Output directory for singularity containers and reference files",
)

OPTION_VERSION = click.option(
    "-v",
    "--container-version",
    show_default=True,
    default=ContainerVersion.RELEASE,
    type=click.Choice(
        [ContainerVersion.DEVELOP, ContainerVersion.RELEASE, balsamic_version]
    ),
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
