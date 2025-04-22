"""Balsamic cluster and submission specific constants."""
from enum import StrEnum
from typing import List

MAX_JOBS: int = 999


class ClusterConfigType(StrEnum):
    """Analysis workflow config type."""

    CACHE: str = "cluster_cache"


class ClusterAccount(StrEnum):
    """Cluster job submission account."""

    DEVELOPMENT: str = "development"


class QOS(StrEnum):
    """Cluster quality of service."""

    LOW: str = "low"
    NORMAL: str = "normal"
    HIGH: str = "high"
    EXPRESS: str = "express"


QOS_OPTIONS: List[QOS] = [qos for qos in QOS]
