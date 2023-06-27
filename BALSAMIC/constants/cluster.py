"""Balsamic cluster and submission specific constants."""
from BALSAMIC.utils.class_types import StrEnum


class ClusterConfigType(StrEnum):
    """Analysis workflow config type."""

    ANALYSIS: str = "cluster_analysis"
    CACHE: str = "cluster_cache"


class ClusterProfile(StrEnum):
    """Profile to submit jobs to the cluster."""

    SLURM: str = "slurm"
    QSUB: str = "qsub"


class QOS(StrEnum):
    """Cluster quality of service."""

    LOW: str = "low"
    NORMAL: str = "normal"
    HIGH: str = "high"
    EXPRESS: str = "express"


class ClusterMailType(StrEnum):
    """Cluster job mail type notification."""

    ALL: str = "ALL"
    BEGIN: str = "BEGIN"
    END: str = "END"
    FAIL: str = "FAIL"
    NONE: str = "NONE"
    REQUEUE: str = "REQUEUE"
    TIME_LIMIT: str = "TIME_LIMIT"
