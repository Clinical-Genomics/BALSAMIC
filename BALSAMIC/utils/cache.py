"""Utility methods for Balsamic init command."""
from typing import Dict

from BALSAMIC import __version__ as balsamic_version
from BALSAMIC.constants.cache import ContainerVersion, DOCKER_URL, DockerContainers


def get_containers(container_version: ContainerVersion) -> Dict[str, str]:
    """Return a dictionary mapping container names to their docker image paths."""
    return {
        container: f"{DOCKER_URL}:{get_docker_image_name(container_version)}-{container}"
        for container in set(DockerContainers)
    }


def get_docker_image_name(container_version: ContainerVersion) -> str:
    """Return docker image base name."""
    if container_version == ContainerVersion.DEVELOP:
        return container_version
    return f"release_v{get_release_version(container_version)}"


def get_release_version(
    container_version: ContainerVersion, current_version: str = balsamic_version
) -> str:
    """Return the containers release version to be downloaded."""
    major, minor, patch = map(int, current_version.split("."))
    if container_version == ContainerVersion.PATCH:
        patch += 1
    elif container_version == ContainerVersion.MINOR:
        minor += 1
        patch = 0
    elif container_version == ContainerVersion.MAJOR:
        major += 1
        minor = 0
        patch = 0
    return f"{major}.{minor}.{patch}"
