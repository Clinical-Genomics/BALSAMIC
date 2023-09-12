"""Helper functions."""
from pathlib import Path
from typing import Dict


def remove_unnecessary_spaces(string: str) -> str:
    """Return a string removing unnecessary empty spaces."""
    return " ".join(string.split())


def get_relative_paths_dict(base_path: Path, data: Dict[str, Path]) -> Dict[str, Path]:
    """Return a dictionary containing relative paths with respect to a given base path."""
    return {key: path.relative_to(base_path) for key, path in data.items()}
