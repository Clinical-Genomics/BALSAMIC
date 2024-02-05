"""Utility methods for Balsamic delivery command."""
import logging
from pathlib import Path
from typing import Dict, Any, List, Generator

from BALSAMIC.utils.exc import BalsamicError

from BALSAMIC.constants.constants import FileType

LOG = logging.getLogger(__name__)


def get_multiqc_deliverables(case_id: str, multiqc_dir: Path) -> List[Dict[str, Any]]:
    """Return a list of MultiQC deliverable files from a directory."""
    multiqc_deliverables: List[Dict[str, Any]] = []
    json_files: Generator[Path, None, None] = multiqc_dir.glob(f"*.{FileType.JSON}")
    for file in json_files:
        deliverable: Dict[str, Any] = {
            "path": file.as_posix(),
            "step": "multiqc",
            "format": FileType.JSON.value,
            "tag": get_file_tags_from_name(file),
            "id": case_id,
        }
        multiqc_deliverables.append(deliverable)
    if not multiqc_deliverables:
        LOG.error(f"No MultiQC deliverable files found in {multiqc_dir.as_posix()}.")
        raise BalsamicError
    return multiqc_deliverables


def get_file_tags_from_name(file_path: Path) -> List[str]:
    """Return Housekeeper tags from the file name after discarding the suffix."""
    return file_path.stem.split("_")
