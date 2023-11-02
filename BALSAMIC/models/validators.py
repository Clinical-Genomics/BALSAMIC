"""Model class validators."""
from pathlib import Path
from typing import Optional


def is_file(file_path: Optional[str]) -> str:
    """Validate file path existence."""
    if file_path:
        if Path(file_path).is_file():
            return file_path
        raise ValueError(f"The supplied file path {file_path} does not exist")
    return file_path


def is_dir(dir_path: Optional[str]) -> str:
    """Validate directory path existence."""
    if dir_path:
        if Path(dir_path).is_dir():
            return dir_path
        raise ValueError(f"The supplied directory path {dir_path} does not exist")
    return dir_path
