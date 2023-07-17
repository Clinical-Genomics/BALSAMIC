"""Balsamic commands common options."""
import click

from BALSAMIC.constants.analysis import RunMode
from BALSAMIC.constants.cache import GenomeVersion
from BALSAMIC.constants.cluster import ClusterProfile, QOS, ClusterMailType

OPTION_GENOME_VERSION = click.option(
    "-g",
    "--genome-version",
    show_default=True,
    default=GenomeVersion.HG19,
    type=click.Choice([GenomeVersion.HG19, GenomeVersion.HG38, GenomeVersion.CanFam3]),
    help="Type and build version of the reference genome",
)

OPTION_RUN_MODE = click.option(
    "--run-mode",
    show_default=True,
    default=RunMode.CLUSTER,
    type=click.Choice([RunMode.CLUSTER, RunMode.LOCAL]),
    help="Run mode to execute Balsamic workflows",
)

OPTION_CLUSTER_PROFILE = click.option(
    "-p",
    "--profile",
    show_default=True,
    default=ClusterProfile.SLURM,
    type=click.Choice([ClusterProfile.SLURM, ClusterProfile.QSUB]),
    help="Cluster profile to submit jobs",
)

OPTION_CLUSTER_QOS = click.option(
    "--qos",
    show_default=True,
    default=QOS.LOW,
    type=click.Choice([QOS.LOW, QOS.NORMAL, QOS.HIGH, QOS.EXPRESS]),
    help="QOS for cluster jobs",
)

OPTION_CLUSTER_ACCOUNT = click.option(
    "--account", type=click.STRING, help="Cluster account to run jobs"
)

OPTION_CLUSTER_MAIL = click.option(
    "--mail-user",
    type=click.STRING,
    help="User email to receive notifications from the cluster",
)

OPTION_CLUSTER_MAIL_TYPE = click.option(
    "--mail-type",
    type=click.Choice(
        [
            ClusterMailType.ALL,
            ClusterMailType.BEGIN,
            ClusterMailType.END,
            ClusterMailType.FAIL,
            ClusterMailType.NONE,
            ClusterMailType.REQUEUE,
            ClusterMailType.TIME_LIMIT,
        ]
    ),
    help="The mail type triggering cluster emails",
)

OPTION_SNAKEMAKE_OPT = click.option(
    "--snakemake-opt",
    multiple=True,
    help="Options to be passed to snakemake",
)

OPTION_FORCE_ALL = click.option(
    "--force-all",
    show_default=True,
    default=False,
    is_flag=True,
    help="Force execution",
)

OPTION_RUN_ANALYSIS = click.option(
    "-r",
    "--run-analysis",
    show_default=True,
    default=False,
    is_flag=True,
    help="Flag to run the actual analysis",
)

OPTION_QUIET = click.option(
    "-q",
    "--quiet",
    default=False,
    is_flag=True,
    help="Instruct snakemake to not output any progress or rule information",
)
