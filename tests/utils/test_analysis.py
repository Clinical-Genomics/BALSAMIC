"""Test Balsamic analysis utility methods."""
from typing import Dict, Any, List

from BALSAMIC.models.cache import CacheConfig
from BALSAMIC.models.snakemake import SingularityBindPath
from BALSAMIC.utils.analysis import (
    get_singularity_bind_paths,
    get_cache_singularity_bind_paths,
)


def test_get_singularity_bind_paths(sample_config: Dict[str, Any]):
    """Test singularity bind paths retrieval."""

    # GIVEN a sample config dictionary

    # WHEN extracting the singularity bind paths
    bind_paths: List[SingularityBindPath] = get_singularity_bind_paths(sample_config)

    # THEN a list of singularity bind paths should be returned
    assert bind_paths
    assert isinstance(bind_paths[0], SingularityBindPath)


def test_get_cache_singularity_bind_paths(cache_config: CacheConfig):
    """Test singularity bind paths retrieval for Balsamic init workflow."""

    # GIVEN a sample config dictionary

    # WHEN extracting the singularity bind paths for balsamic init workflow
    bind_paths: List[SingularityBindPath] = get_cache_singularity_bind_paths(
        cache_config
    )

    # THEN a list of cache singularity bind paths should be returned
    assert bind_paths
    assert isinstance(bind_paths[0], SingularityBindPath)
