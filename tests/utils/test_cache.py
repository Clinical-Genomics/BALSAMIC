"""Test init utility methods."""
from typing import Dict

import pytest

from BALSAMIC import __version__ as current_version
from BALSAMIC.constants.cache import ContainerVersion
from BALSAMIC.utils.cache import (
    get_containers,
    get_docker_image_name,
    get_release_version,
)


def test_get_containers(develop_containers: Dict[str, str]):
    """Test containers retrieval given a container version."""

    # GIVEN a container version

    # WHEN getting the containers dictionary
    containers: Dict[str, str] = get_containers(ContainerVersion.DEVELOP)

    # THEN the version associated containers should be returned
    assert containers == develop_containers


@pytest.mark.parametrize(
    "container_version, expected_image_name",
    [
        (ContainerVersion.DEVELOP, "develop"),
        (ContainerVersion.RELEASE, f"release_v{current_version}"),
    ],
)
def test_get_docker_image_name(
    container_version: ContainerVersion, expected_image_name: str
):
    """Test container name retrieval given a container version."""

    # GIVEN a container version

    # WHEN getting the image name
    image_name: str = get_docker_image_name(container_version)

    # THEN the correct image name should be returned
    assert image_name == expected_image_name


@pytest.mark.parametrize(
    "container_version, expected_release_version",
    [
        (ContainerVersion.PATCH, "1.0.1"),
        (ContainerVersion.MINOR, "1.1.0"),
        (ContainerVersion.MAJOR, "2.0.0"),
        (ContainerVersion.RELEASE, "1.0.0"),
    ],
)
def test_get_release_version(
    balsamic_version: str,
    container_version: ContainerVersion,
    expected_release_version: str,
):
    """Test release version retrieval."""

    # GIVEN a mocked balsamic version

    # WHEN retrieving the containers version
    release_version: str = get_release_version(container_version, balsamic_version)

    # THEN the correct version should be returned
    assert release_version == expected_release_version
