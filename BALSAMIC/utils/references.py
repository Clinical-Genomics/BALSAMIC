"""Utility functions for Balsamic references."""
from pathlib import Path


def merge_reference_metadata(
    existing_refs: dict, observation_paths: dict, metadata_lookup: dict
) -> dict:
    """
    Merge user-provided observation file paths and existing references
    with static metadata definitions.
    Returns a new merged references dict.
    """
    merged = {}

    # 1. Start with existing references
    for key, value in existing_refs.items():
        if isinstance(value, (str, Path)):
            entry = {"file": Path(value)}
        elif isinstance(value, dict):
            entry = {**value, "file": Path(value["file"])}
        else:
            raise TypeError(f"Unsupported reference format for '{key}': {value!r}")

        if key in metadata_lookup:
            # Merge in static metadata without overwriting file path
            for meta_key, meta_val in metadata_lookup[key].items():
                entry.setdefault(meta_key, meta_val)

        merged[key] = entry

    # 2. Add/update with command-line–provided observation paths
    for key, file_path in observation_paths.items():
        if file_path is None:
            continue

        entry = {"file": Path(file_path)}
        if key in metadata_lookup:
            entry.update(metadata_lookup[key])

        merged[key] = entry  # overrides existing if same key

    return merged
