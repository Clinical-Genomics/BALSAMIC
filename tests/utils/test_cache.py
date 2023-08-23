"""Test init utility methods."""
from typing import Dict

from BALSAMIC.utils.cache import get_containers


def test_get_containers_develop(develop_containers: Dict[str, str]):
    """Test containers retrieval given a develop cache version."""

    # GIVEN a cache version

    # WHEN getting the containers dictionary
    containers: Dict[str, str] = get_containers("develop")

    # THEN the version associated containers should be returned
    assert containers == develop_containers
