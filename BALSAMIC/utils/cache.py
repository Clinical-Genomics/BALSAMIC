"""Utility methods for Balsamic init command."""
from typing import Dict

from BALSAMIC.constants.cache import DOCKER_URL, DockerContainers, CacheVersion


def get_containers(cache_version: str) -> Dict[str, str]:
    """Return a dictionary mapping container names to their docker image paths."""
    cache_version: str = (
        cache_version
        if cache_version == CacheVersion.DEVELOP
        else f"release_v{cache_version}"
    )
    return {
        container: f"{DOCKER_URL}:{cache_version}-{container}"
        for container in set(DockerContainers)
    }
