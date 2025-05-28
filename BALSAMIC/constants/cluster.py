"""Balsamic cluster and submission specific constants."""
from enum import StrEnum
from typing import List

MAX_JOBS: int = 999


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
