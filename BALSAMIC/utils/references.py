"""Utility functions for Balsamic references."""
from pathlib import Path
from BALSAMIC.constants.analysis import VARIANT_OBSERVATION_METAVALUES
from typing import Dict, Any


def merge_reference_metadata(
    existing_refs: dict[str, str],
    observation_paths: dict[str, str] | None = None,
) -> dict[str, dict]:
    """
    Build a merged references dict from simple {key: filepath} inputs.

    - existing_refs: baseline paths; each becomes {"file": Path(...)} plus static metadata if available.
    - observation_paths: optional overrides/additions applied after existing_refs; None values are ignored.

    Returns:
        { key: {"file": Path(...), **static_meta_if_any} }
    """
    merged: Dict[str, Dict[str, Any]] = {}

    # 1) Seed with existing references
    for key, fp in existing_refs.items():
        entry: Dict[str, Any] = {"file": Path(fp)}
        entry.update(VARIANT_OBSERVATION_METAVALUES.get(key, {}))
        merged[key] = entry

    # 2) Apply observation overrides/additions
    if observation_paths:
        for key, fp in observation_paths.items():
            if fp is None:
                continue
            entry: Dict[str, Any] = {"file": Path(fp)}
            entry.update(VARIANT_OBSERVATION_METAVALUES.get(key, {}))
            merged[key] = entry  # override or add

    return merged
