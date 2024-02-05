"""Test Balsamic delivery utility methods."""
from pathlib import Path
from typing import Any, Dict, List

import pytest
from BALSAMIC.utils.exc import BalsamicError

from BALSAMIC.constants.constants import FileType
from BALSAMIC.utils.delivery import get_file_tags_from_name, get_multiqc_deliverables


def test_get_multiqc_deliverables(case_id_tumor_only: str, multiqc_data_dir: Path):
    """Test MultiQC delivery files parsing."""

    # GIVEN a case ID and a MultiQC data directory

    # WHEN extracting the deliverables
    multiqc_deliverables: List[Dict[str, Any]] = get_multiqc_deliverables(
        case_id=case_id_tumor_only, multiqc_dir=multiqc_data_dir
    )

    # THEN the correct number of deliverables should be returned with the expected structure
    assert len(multiqc_deliverables) == 5
    assert all(isinstance(item["path"], str) for item in multiqc_deliverables)
    assert all(item["step"] == "multiqc" for item in multiqc_deliverables)
    assert all(item["format"] == FileType.JSON for item in multiqc_deliverables)
    assert all(isinstance(item["tag"], list) for item in multiqc_deliverables)
    assert all(item["id"] == case_id_tumor_only for item in multiqc_deliverables)


def test_get_multiqc_deliverables_error(
    case_id_tumor_only: str, fastq_dir_tumor_only: str
):
    """Test MultiQC delivery files parsing when incorrect path is provided."""

    # GIVEN a case ID and an incorrect MultiQC data directory

    # WHEN extracting the deliverables
    with pytest.raises(BalsamicError):
        # THEN an exception should be raised
        get_multiqc_deliverables(
            case_id=case_id_tumor_only, multiqc_dir=Path(fastq_dir_tumor_only)
        )


def test_get_file_tags_from_name():
    """Test tag extraction from a file name."""

    # GIVEN a mock file object
    file_path: Path = Path(
        "/analysis/qc/multiqc_data/multiqc_picard_AlignmentSummaryMetrics.json"
    )

    # WHEN extracting the tags from the file name
    tags: List[str] = get_file_tags_from_name(file_path)

    # THEN the correct tags are extracted from the file name
    assert tags == ["multiqc", "picard", "AlignmentSummaryMetrics"]
