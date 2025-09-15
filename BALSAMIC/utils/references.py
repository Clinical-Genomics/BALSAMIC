"""Utility functions for Balsamic references."""
from pathlib import Path
from BALSAMIC.constants.analysis import VARIANT_OBSERVATION_METAVALUES
from typing import Dict, Any


def add_reference_metadata(
    references: dict[str, str],
) -> dict[str, dict]:
    """
    Build a rich reference map by attaching file Paths and static metadata.

    For each input entry {key: filepath}, returns:
        { key: {"file": Path(filepath), **VARIANT_OBSERVATION_METAVALUES.get(key, {})} }

    The input dict is not modified.
    """
    references_metadata: Dict[str, Dict[str, Any]] = {}

    for key, fp in references.items():
        entry: Dict[str, Any] = {"file": Path(fp)}
        entry.update(VARIANT_OBSERVATION_METAVALUES.get(key, {}))
        references_metadata[key] = entry

    return references_metadata
