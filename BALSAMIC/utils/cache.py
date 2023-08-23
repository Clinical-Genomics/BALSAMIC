"""Utility methods for Balsamic init command."""
from typing import Dict

from BALSAMIC.constants.cache import DOCKER_URL, DockerContainers


def get_containers(cache_version: str) -> Dict[str, str]:
    """Return a dictionary mapping container names to their docker image paths."""
    return {
        container: f"{DOCKER_URL}:{cache_version}-{container}"
        for container in set(DockerContainers)
    }
