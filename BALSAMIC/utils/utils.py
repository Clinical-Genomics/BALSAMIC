"""Helper functions."""
from pathlib import Path
from typing import Dict


def remove_unnecessary_spaces(string: str) -> str:
    """Return a string removing unnecessary empty spaces."""
    return " ".join(string.split())


def get_relative_paths_dict(base_path: Path, data: Dict[str, Path]) -> Dict[str, str]:
    """Return a dictionary containing relative paths with respect to a given base path."""
    return {key: path.relative_to(base_path).as_posix() for key, path in data.items()}


def get_absolute_paths_dict(base_path: Path, data: Dict[str, Path]) -> Dict[str, Path]:
    """Return a dictionary containing absolute resolved paths with respect to a given base path."""
    return {key: Path(base_path, path).resolve() for key, path in data.items()}
