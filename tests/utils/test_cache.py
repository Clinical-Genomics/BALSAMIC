"""Test init utility methods."""
from typing import Dict

from BALSAMIC.constants.cache import ContainerVersion

from BALSAMIC.utils.cache import get_containers


def test_get_containers(develop_containers: Dict[str, str]):
    """Test containers retrieval given a container version."""

    # GIVEN a container version

    # WHEN getting the containers dictionary
    containers: Dict[str, str] = get_containers(ContainerVersion.DEVELOP)

    # THEN the version associated containers should be returned
    assert containers == develop_containers
