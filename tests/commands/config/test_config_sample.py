import os
import json
import pytest

from BALSAMIC.commands.config import get_config, write_json
from BALSAMIC.tools import iterdict


def test_get_config(config_files):
    # GIVEN the config files name
    files = [
        "install", "reference", "sample", "analysis_paired",
        "analysis_paired_umi", "analysis_single", "analysis_single_umi"
    ]
    # WHEN passing file names
    for file in files:
        #THEN return the config files path
        assert config_files[file] in get_config(file)


def test_write_json(tmp_path, config_files):
    # GIVEN a dict from sample json file (reference.json)
    ref_json = json.load(open(config_files['reference'], 'r'))

    tmp = tmp_path / "tmp"
    tmp.mkdir()
    output_json = tmp / "output.json"

    #WHEN passing dict and file name
    write_json(ref_json, output_json)
    output = output_json.read_text()

    #THEN It will create a json file with given dict
    for key, value in iterdict(ref_json):
        assert key in output
        assert value in output

    assert len(list(tmp.iterdir())) == 1
