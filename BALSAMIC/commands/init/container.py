import re
import logging
import subprocess
from pathlib import Path
import click

from BALSAMIC import __version__ as balsamic_version
from BALSAMIC.utils.exc import BalsamicError

LOG = logging.getLogger(__name__)


@click.command("container",
               short_help="Download matching version for container")
@click.option("-o",
              "--out-dir",
              required=True,
              help="Output directory for container files.")
@click.option("-v",
              "--container-version",
              show_default=True,
              default=balsamic_version,
              help="Container for BALSAMIC version to download")
@click.option('-f',
              '--force',
              show_default=True,
              default=False,
              is_flag=True,
              help="Force re-downloading all containers")
@click.pass_context
def container(context, container_version, out_dir, force):
    """
    Pull container(s) for BALSAMIC according to matching version
    """
    # resolve out_dir to absolute path
    out_dir = Path(out_dir).resolve()

    pattern = re.compile(r"^(\d+\.)?(\d+\.)?(\*|\d+)$")
    if pattern.findall(container_version):
        docker_image_name = "release_v" + container_version
    else:
        docker_image_name = container_version

    container_stub_url = "docker://hassanf/balsamic:" + docker_image_name

    # Pull container
    LOG.info("Pulling singularity image {}.".format(container_stub_url))

    # Set container name
    image_name = Path(out_dir,
                      "BALSAMIC_{}.sif".format(docker_image_name)).as_posix()
    LOG.info("Image will be downloaded to {}".format(image_name))

    try:
        subprocess.check_output(
            [
                "singularity",
                "pull",
                "--name",
                "{}".format(image_name),
                container_stub_url,
            ],
            cwd=out_dir,
            stderr=subprocess.STDOUT,
        )
    except subprocess.CalledProcessError as e:
        raise BalsamicError("Failed to pull singularity image "
                            "from {}:\n{}".format(container_stub_url,
                                                  e.stdout.decode()))
