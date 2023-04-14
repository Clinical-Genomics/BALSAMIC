"""Balsamic cache specific constants."""
from BALSAMIC.utils.str_enum import StrEnum


class ContainerVersion(StrEnum):
    """Balsamic container versions."""

    DEVELOP: str = "develop"
    RELEASE: str = "release"


class FileType(StrEnum):
    """Balsamic reference file types."""

    BED: str = "bed"
    DICT: str = "dict"
    FAI: str = "fai"
    FASTA: str = "fasta"
    FLAT: str = "flat"
    GFF: str = "gff"
    GTF: str = "gtf"
    GZ: str = "gz"
    TEXT: str = "text"
    TSV: str = "tsv"
    TXT: str = "txt"
    SIF: str = "sif"
    VCF: str = "vcf"


class BwaIndexFileType(StrEnum):
    """BWA genome index file suffixes."""

    AMB: str = "amb"
    ANN: str = "ann"
    BWT: str = "bwt"
    PAC: str = "pac"
    SA: str = "sa"


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
