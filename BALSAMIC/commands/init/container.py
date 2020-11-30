import re
import logging
import subprocess
from pathlib import Path
import click

from BALSAMIC import __version__ as balsamic_version
from BALSAMIC.utils.constants import BALSAMIC_DOCKER_PATH
from BALSAMIC.utils.exc import BalsamicError

LOG = logging.getLogger(__name__)


@click.command("container",
               short_help="Download matching version for container")
@click.option("-d",
              "--dry",
              show_default=True,
              default=False,
              is_flag=True,
              help="Dry run mode.")
@click.option("-o",
              "--out-dir",
              required=True,
              help="Output directory for container files.")
@click.option("-v",
              "--container-version",
              show_default=True,
              default=balsamic_version,
              type=click.Choice(["develop", "master", balsamic_version]),
              help="Container for BALSAMIC version to download")
@click.option('-f',
              '--force',
              show_default=True,
              default=False,
              is_flag=True,
              help="Force re-downloading all containers")
@click.pass_context
def container(context, container_version, out_dir, force, dry):
    """
    Pull container(s) for BALSAMIC according to matching version
    """
    LOG.info("BALSAMIC started with log level %s" % context.obj['loglevel'])
    # resolve out_dir to absolute path
    out_dir = Path(out_dir).resolve()

    pattern = re.compile(r"^(\d+\.)?(\d+\.)?(\*|\d+)$")
    if pattern.findall(container_version):
        docker_image_name = "release_v{}".format(container_version)
    else:
        docker_image_name = container_version

    container_stub_url = "{}:{}".format(BALSAMIC_DOCKER_PATH,
                                        docker_image_name)

    # Pull container
    LOG.info("Singularity image source: {}".format(container_stub_url))

    # Set container name according to above docker image name
    image_name = Path(out_dir,
                      "BALSAMIC_{}.sif".format(docker_image_name)).as_posix()
    LOG.info("Image will be downloaded to {}".format(image_name))
    LOG.info("Starting download. This process can take some time...")

    cmd = ["singularity", "pull", "--name", f"{image_name}"]
    if force:
        cmd.append("--force")
    cmd.append(container_stub_url)

    try:
        if dry:
            LOG.info("Dry run mode, The following command will run: {}".format(
                " ".join(cmd)))
        else:
            subprocess.check_output(cmd, cwd=out_dir, stderr=subprocess.STDOUT)

    except:
        LOG.error("Failed to pull singularity image "
                  "from {}".format(container_stub_url))
        raise click.Abort()
