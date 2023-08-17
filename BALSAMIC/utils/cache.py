"""Utility methods for Balsamic init command."""
from typing import Dict

from BALSAMIC import __version__ as balsamic_version
from BALSAMIC.constants.cache import CacheVersion, DOCKER_URL, DockerContainers


def get_containers(cache_version: CacheVersion) -> Dict[str, str]:
    """Return a dictionary mapping container names to their docker image paths."""
    return {
        container: f"{DOCKER_URL}:{get_docker_image_name(cache_version)}-{container}"
        for container in set(DockerContainers)
    }


def get_docker_image_name(cache_version: CacheVersion) -> str:
    """Return docker image base name."""
    if cache_version == CacheVersion.DEVELOP:
        return cache_version
    return f"release_v{get_release_version(cache_version)}"


def get_release_version(
    cache_version: CacheVersion, current_version: str = balsamic_version
) -> str:
    """Return the containers release version to be downloaded."""
    major, minor, patch = map(int, current_version.split("."))
    if cache_version == CacheVersion.PATCH:
        patch += 1
    elif cache_version == CacheVersion.MINOR:
        minor += 1
        patch = 0
    elif cache_version == CacheVersion.MAJOR:
        major += 1
        minor = 0
        patch = 0
    return f"{major}.{minor}.{patch}"
