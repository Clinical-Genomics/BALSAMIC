"""Test init utility methods."""
from typing import Dict

import pytest

from BALSAMIC import __version__ as current_version
from BALSAMIC.constants.cache import CacheVersion
from BALSAMIC.utils.cache import (
    get_containers,
    get_docker_image_name,
    get_release_version,
)


def test_get_containers(develop_containers: Dict[str, str]):
    """Test containers retrieval given a cache version."""

    # GIVEN a cache version

    # WHEN getting the containers dictionary
    containers: Dict[str, str] = get_containers(CacheVersion.DEVELOP)

    # THEN the version associated containers should be returned
    assert containers == develop_containers


@pytest.mark.parametrize(
    "cache_version, expected_image_name",
    [
        (CacheVersion.DEVELOP, "develop"),
        (CacheVersion.RELEASE, f"release_v{current_version}"),
    ],
)
def test_get_docker_image_name(cache_version: CacheVersion, expected_image_name: str):
    """Test container name retrieval given a cache version."""

    # GIVEN a cache version

    # WHEN getting the image name
    image_name: str = get_docker_image_name(cache_version)

    # THEN the correct image name should be returned
    assert image_name == expected_image_name


@pytest.mark.parametrize(
    "cache_version, expected_release_version",
    [
        (CacheVersion.PATCH, "1.0.1"),
        (CacheVersion.MINOR, "1.1.0"),
        (CacheVersion.MAJOR, "2.0.0"),
        (CacheVersion.RELEASE, "1.0.0"),
    ],
)
def test_get_release_version(
    balsamic_version: str,
    cache_version: CacheVersion,
    expected_release_version: str,
):
    """Test release version retrieval."""

    # GIVEN a mocked balsamic version

    # WHEN retrieving the cache version
    release_version: str = get_release_version(cache_version, balsamic_version)

    # THEN the correct version should be returned
    assert release_version == expected_release_version
