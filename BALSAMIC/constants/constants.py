"""General use constants."""
from BALSAMIC.utils.class_types import StrEnum

EXIT_SUCCESS: int = 0
EXIT_FAIL: int = 1


class FileType(StrEnum):
    """Balsamic analysis and reference file extensions."""

    BED: str = "bed"
    DICT: str = "dict"
    FAI: str = "fai"
    FASTA: str = "fasta"
    FLAT: str = "flat"
    GFF: str = "gff"
    GTF: str = "gtf"
    GZ: str = "gz"
    JSON: str = "json"
    LOG: str = "log"
    SIF: str = "sif"
    TBI: str = "tbi"
    TSV: str = "tsv"
    TXT: str = "txt"
    VCF: str = "vcf"


class BwaIndexFileType(StrEnum):
    """BWA genome index file extensions."""

    AMB: str = "amb"
    ANN: str = "ann"
    BWT: str = "bwt"
    PAC: str = "pac"
    SA: str = "sa"
