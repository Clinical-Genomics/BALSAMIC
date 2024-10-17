"""Script to submit jobs to a cluster."""
import shutil
from typing import Any, Dict, List, Optional

import click
from snakemake import utils

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
@click.argument("dependencies", nargs=-1, type=click.STRING)
@click.argument("job_script", nargs=1, type=click.Path(exists=True, resolve_path=True))
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
    case_id: str,
    job_script: str,
    log_dir: str,
    profile: ClusterProfile,
    script_dir: str,
    benchmark: Optional[bool] = False,
    dependencies: Optional[List[str]] = None,
    mail_type: Optional[str] = None,
    mail_user: Optional[str] = None,
    qos: Optional[QOS] = QOS.LOW,
) -> None:
    """
    Submits jobs to the cluster. Each job is submitted sequentially, and their respective job IDs are collected
    from the output. These job IDs are then forwarded as dependencies to the subsequent jobs.
    """
    job_script: str = shutil.copy2(src=job_script, dst=script_dir)
    job_properties: Dict[str, Any] = utils.read_job_properties(job_script)
    scheduler: Scheduler = Scheduler(
        account=account,
        benchmark=benchmark,
        case_id=case_id,
        dependencies=dependencies,
        job_properties=job_properties,
        job_script=job_script,
        log_dir=log_dir,
        mail_type=mail_type,
        mail_user=mail_user,
        profile=profile,
        qos=qos,
    )
    scheduler.submit_job()


if __name__ == "__main__":
    immediate_submit()
