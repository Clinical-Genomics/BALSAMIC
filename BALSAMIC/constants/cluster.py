"""Balsamic cluster and submission specific constants."""
from BALSAMIC.utils.class_types import StrEnum


class ClusterConfigType(StrEnum):
    """Analysis workflow config type."""

    ANALYSIS: str = "cluster_analysis"
    CACHE: str = "cluster_cache"
