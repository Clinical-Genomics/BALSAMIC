"""Script to submit jobs to a cluster."""
import shutil
from typing import Any, Dict, List

import click
from snakemake.utils import read_job_properties

from BALSAMIC.commands.options import (
    OPTION_BENCHMARK,
    OPTION_CLUSTER_ACCOUNT,
    OPTION_CLUSTER_MAIL,
    OPTION_CLUSTER_MAIL_TYPE,
    OPTION_CLUSTER_PROFILE,
    OPTION_CLUSTER_QOS,
)
from BALSAMIC.constants.cluster import QOS, ClusterProfile
from BALSAMIC.models.scheduler import Scheduler


@click.command()
@click.argument("case_id", nargs=1, required=True, type=click.STRING)
@click.argument(
    "dependencies",
    nargs=-1,
    type=click.STRING,
    help="Snakemake job dependencies",
)
@click.argument(
    "job_script",
    nargs=1,
    type=click.Path(exists=True, resolve_path=True),
    help="Snakemake job script path",
)
@OPTION_CLUSTER_ACCOUNT
@OPTION_BENCHMARK
@OPTION_CLUSTER_MAIL_TYPE
@OPTION_CLUSTER_MAIL
@OPTION_CLUSTER_PROFILE
@OPTION_CLUSTER_QOS
@click.option(
    "--log-dir",
    type=click.Path(exists=True, resolve_path=True),
    required=True,
    help="Logging directory path",
)
@click.option(
    "--script-dir",
    type=click.Path(exists=True, resolve_path=True),
    required=True,
    help="Script directory path",
)
def immediate_submit(
    account: str,
    benchmark: bool,
    case_id: str,
    dependencies: List[str],
    job_script: str,
    log_dir: str,
    mail_type: str,
    mail_user: str,
    profile: ClusterProfile,
    qos: QOS,
    script_dir: str,
) -> None:
    """
    Submits jobs to the cluster. Each job is submitted sequentially, and their respective job IDs are collected
    from the output. These job IDs are then forwarded as dependencies to the subsequent jobs.
    """
    job_script: str = shutil.copy2(src=job_script, dst=script_dir)
    job_properties: Dict[str, Any] = read_job_properties(job_script)
    scheduler: Scheduler = Scheduler(
        account=account,
        benchmark=benchmark,
        dependencies=dependencies,
        job_properties=job_properties,
        job_script=job_script,
        log_dir=log_dir,
        mail_type=mail_type,
        mail_user=mail_user,
        profile=profile,
        qos=qos,
    )


if __name__ == "__main__":
    immediate_submit()
