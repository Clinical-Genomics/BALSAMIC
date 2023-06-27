"""Balsamic init workflow fixtures."""
from typing import Dict

import pytest

from BALSAMIC.constants.cache import DockerContainers


@pytest.fixture(scope="session", name="develop_containers")
def fixture_develop_containers() -> Dict[str, str]:
    """Develop containers fixture."""
    return {
        DockerContainers.ASCAT.value: "docker://clinicalgenomics/balsamic:develop-ascatNgs",
        DockerContainers.VCF2CYTOSURE.value: "docker://clinicalgenomics/balsamic:develop-vcf2cytosure",
        DockerContainers.PYTHON_3.value: "docker://clinicalgenomics/balsamic:develop-varcall_py3",
        DockerContainers.BALSAMIC.value: "docker://clinicalgenomics/balsamic:develop-balsamic",
        DockerContainers.SOMALIER.value: "docker://clinicalgenomics/balsamic:develop-somalier",
        DockerContainers.CNVPYTOR.value: "docker://clinicalgenomics/balsamic:develop-cnvpytor",
        DockerContainers.ALIGN_QC.value: "docker://clinicalgenomics/balsamic:develop-align_qc",
        DockerContainers.ANNOTATE.value: "docker://clinicalgenomics/balsamic:develop-annotate",
        DockerContainers.PYTHON_27.value: "docker://clinicalgenomics/balsamic:develop-varcall_py27",
        DockerContainers.CNVKIT.value: "docker://clinicalgenomics/balsamic:develop-varcall_cnvkit",
        DockerContainers.COVERAGE_QC.value: "docker://clinicalgenomics/balsamic:develop-coverage_qc",
        DockerContainers.DELLY.value: "docker://clinicalgenomics/balsamic:develop-delly",
    }
