"""Input/Output utils file"""

import json
from pathlib import Path

import yaml


def read_json(json_path) -> dict:
    if Path(json_path).exists():
        with open(json_path, "r") as fn:
            return json.load(fn)
    else:
        raise FileNotFoundError(f"The JSON file {json_path} was not found.")


def write_json(json_out, output_config):
    """Writes JSON format data to an output file"""
    try:
        with open(output_config, "w") as fn:
            json.dump(json_out, fn, indent=4)
    except OSError as error:
        raise OSError(f"Error while writing JSON file: {output_config}, error: {error}")


def read_yaml(yaml_path):
    """Retrieves data from a yaml file"""
    if Path(yaml_path).exists():
        with open(yaml_path, "r") as fn:
            return yaml.load(fn, Loader=yaml.SafeLoader)
    else:
        raise FileNotFoundError(f"The YAML file {yaml_path} was not found.")
