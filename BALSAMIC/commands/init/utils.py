"""Utility methods for Balsamic init command."""
from typing import Dict

from BALSAMIC import __version__ as balsamic_version

from BALSAMIC.constants.cache import ContainerVersion, DOCKER_URL, DockerContainers


def get_containers(container_version: ContainerVersion) -> Dict[str, str]:
    """Return a dictionary mapping container names to their docker image paths."""
    image_name: str = (
        f"release_v{balsamic_version}"
        if container_version == ContainerVersion.RELEASE
        else container_version
    )
    return {
        container: f"{DOCKER_URL}:{image_name}-{container}"
        for container in set(DockerContainers)
    }
