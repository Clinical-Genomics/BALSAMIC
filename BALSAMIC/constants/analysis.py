from BALSAMIC.constants.common import StrEnum


class GenomeVersion(StrEnum):
    """Reference genome version"""

    HG19: str = "hg19"
    HG38: str = "hg38"
    CanFam3: str = "canfam3"


class RunMode(StrEnum):
    """Run mode to generate the balsamic cache."""

    CLUSTER: str = "cluster"
    LOCAL: str = "local"


class ClusterProfile(StrEnum):
    """Profile to submit jobs to the cluster."""

    SLURM: str = "slurm"
    QSUB: str = "qsub"


class QOS(StrEnum):
    """Cluster QOS."""

    LOW: str = "low"
    NORMAL: str = "normal"
    HIGH: str = "high"
    EXPRESS: str = "express"


class ClusterMailType(StrEnum):
    """Cluster mail type."""

    ALL: str = "ALL"
    BEGIN: str = "BEGIN"
    END: str = "END"
    FAIL: str = "FAIL"
    NONE: str = "NONE"
    REQUEUE: str = "REQUEUE"
    TIME_LIMIT: str = "TIME_LIMIT"
