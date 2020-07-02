import os
import re
import json
import pytest
from unittest import mock
import click

from pathlib import Path



def test_dag_graph_success(tumor_normal_wgs_config, tumor_only_config, tumor_normal_config, tumor_only_wgs_config):
    # WHEN creating config using standard CLI input
    # THEN DAG graph should be created successfully
    assert Path(json.load(open(tumor_normal_wgs_config))["analysis"]["dag"]).exists()
    assert Path(json.load(open(tumor_normal_config))["analysis"]["dag"]).exists()
    assert Path(json.load(open(tumor_only_wgs_config))["analysis"]["dag"]).exists()
    assert Path(json.load(open(tumor_only_config))["analysis"]["dag"]).exists()




