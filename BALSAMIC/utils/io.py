"""Input/Output utils file."""

import json
import logging
import shutil
from pathlib import Path

import yaml

LOG = logging.getLogger(__name__)


def read_json(json_path: str) -> dict:
    """Read JSON file and return a dictionary."""
    if Path(json_path).exists():
        with open(json_path, "r") as fn:
            return json.load(fn)
    else:
        raise FileNotFoundError(f"The JSON file {json_path} was not found")


def write_json(json_obj: dict, path: str) -> None:
    """Write JSON format data to an output file."""
    try:
        with open(path, "w") as fn:
            json.dump(json_obj, fn, indent=4)
    except OSError as error:
        raise OSError(f"Error while writing JSON file: {path}, error: {error}")


def read_yaml(yaml_path: str) -> dict:
    """Read data from a yaml file."""
    if Path(yaml_path).exists():
        with open(yaml_path, "r") as fn:
            return yaml.load(fn, Loader=yaml.SafeLoader)
    else:
        raise FileNotFoundError(f"The YAML file {yaml_path} was not found")
