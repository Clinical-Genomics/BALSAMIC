"""Balsamic cluster and submission specific constants."""
from typing import List

from BALSAMIC.utils.class_types import StrEnum


MAX_JOBS: int = 999


class ClusterConfigType(StrEnum):
    """Analysis workflow config type."""

    ANALYSIS: str = "cluster_analysis"
    CACHE: str = "cluster_cache"


class ClusterProfile(StrEnum):
    """Profile to submit jobs to the cluster."""

    SLURM: str = "slurm"
    QSUB: str = "qsub"


CLUSTER_PROFILES: List[ClusterProfile] = [profile.value for profile in ClusterProfile]


class ClusterAccount(StrEnum):
    """Cluster job submission account."""

    DEVELOPMENT: str = "development"


class QOS(StrEnum):
    """Cluster quality of service."""

    LOW: str = "low"
    NORMAL: str = "normal"
    HIGH: str = "high"
    EXPRESS: str = "express"


QOS_OPTIONS: List[QOS] = [qos.value for qos in QOS]


class ClusterMailType(StrEnum):
    """Cluster job mail type notification."""

    ALL: str = "ALL"
    BEGIN: str = "BEGIN"
    END: str = "END"
    FAIL: str = "FAIL"
    NONE: str = "NONE"
    REQUEUE: str = "REQUEUE"
    TIME_LIMIT: str = "TIME_LIMIT"


CLUSTER_MAIL_TYPES: List[ClusterMailType] = [type.value for type in ClusterMailType]
