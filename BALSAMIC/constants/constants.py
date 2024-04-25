"""General use constants."""
from enum import StrEnum
from typing import List

EXIT_SUCCESS: int = 0
EXIT_FAIL: int = 1


class LogLevel(StrEnum):
    NOTSET = "NOTSET"
    DEBUG = "DEBUG"
    INFO = "INFO"
    WARNING = "WARNING"
    ERROR = "ERROR"
    FATAL = "FATAL"
    CRITICAL = "CRITICAL"


LOG_LEVELS: List[LogLevel] = [level for level in LogLevel]


class FileType(StrEnum):
    """Balsamic analysis and reference file extensions."""

    BED: str = "bed"
    CSV: str = "csv"
    DICT: str = "dict"
    FAI: str = "fai"
    FASTA: str = "fasta"
    FASTQ: str = "fastq"
    FLAT: str = "flat"
    GFF: str = "gff"
    GTF: str = "gtf"
    GZ: str = "gz"
    JSON: str = "json"
    LOG: str = "log"
    PDF: str = "pdf"
    SIF: str = "sif"
    TBI: str = "tbi"
    TSV: str = "tsv"
    TXT: str = "txt"
    VCF: str = "vcf"
    YAML: str = "yaml"


class BwaIndexFileType(StrEnum):
    """BWA genome index file extensions."""

    AMB: str = "amb"
    ANN: str = "ann"
    BWT: str = "bwt"
    PAC: str = "pac"
    SA: str = "sa"
