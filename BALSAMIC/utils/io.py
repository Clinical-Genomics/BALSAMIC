import json
from pathlib import Path


def read_json(json_path) -> dict:
    if Path(json_path).exists():
        with open(json_path, "r") as fn:
            return json.load(fn)
    else:
        raise FileNotFoundError(f"The JSON file {json_path} was not found.")
