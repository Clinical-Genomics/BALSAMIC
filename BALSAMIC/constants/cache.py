from BALSAMIC.constants.common import StrEnum


class ContainerVersion(StrEnum):
    """Balsamic container versions."""

    DEVELOP: str = "develop"
    RELEASE: str = "release"


DOCKER_PATH: str = "docker://clinicalgenomics/balsamic"

DOCKER_CONTAINERS: set = {
    "align_qc",
    "annotate",
    "ascatNgs",
    "balsamic",
    "cnvpytor",
    "coverage_qc",
    "delly",
    "somalier",
    "varcall_cnvkit",
    "varcall_py3",
    "varcall_py27",
    "vcf2cytosure",
}
