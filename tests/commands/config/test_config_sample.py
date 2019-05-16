import os
import json
import pytest

from BALSAMIC.commands.config import get_config, write_json


def iterdict(dic):
    """ dictionary iteration - returns generator"""
    for key, value in dic.items():
        if isinstance(value, dict):
            yield from iterdict(value)
        else:
            yield key, value


def test_get_config(config_files):
    """ testing get config function"""
    files = [
        "install", "reference", "sample", "analysis_paired",
        "analysis_paired_umi", "analysis_single", "analysis_single_umi"
    ]

    for file in files:
        assert config_files[file] in get_config(file)


def test_write_json(tmp_path, config_files):
    """ tested by tmp_path fixture using reference.json data """
    ref_json = json.load(open(config_files['reference'], 'r'))

    tmp = tmp_path / "tmp"
    tmp.mkdir()
    output_json = tmp / "output.json"
    write_json(ref_json, output_json)
    output = output_json.read_text()

    for key, value in iterdict(ref_json):
        print(key)
        assert key in output
        assert value in output
