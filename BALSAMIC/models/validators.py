"""Model common validators."""
import logging
from pathlib import Path

from BALSAMIC.constants.analysis import GenomeVersion

LOG = logging.getLogger(__name__)


def validate_dir_path(directory: str) -> str:
    """Check that the provided directory path exists."""
    if not Path(directory).is_dir():
        LOG.error(f"The provided directory path ({directory}) does not exist")
        raise ValueError

    return directory


def validate_file_path(file: str) -> str:
    """Check that the provided directory path exists."""
    if not Path(file).is_file():
        LOG.error(f"The provided file path ({file}) does not exist")
        raise ValueError

    return file


def validate_genome_version(genome_version: GenomeVersion) -> str:
    """Check if the provided genome version is supported by Balsamic."""
    if genome_version not in list(GenomeVersion):
        LOG.error(f"The provided genome version ({genome_version}) is not accepted")
        raise ValueError

    return genome_version
