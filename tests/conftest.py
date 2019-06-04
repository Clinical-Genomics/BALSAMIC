import pytest

from functools import partial
from click.testing import CliRunner

from BALSAMIC.commands import cli


@pytest.fixture
def cli_runner():
    """ click - cli testing """
    runner = CliRunner()
    return runner


@pytest.fixture
def invoke_cli(cli_runner):
    """ invoking cli commands with options"""
    return partial(cli_runner.invoke, cli)


@pytest.fixture(scope='session')
def config_files():
    """ dict: path of the config files """
    return {
        "install": "BALSAMIC/config/install.json",
        "reference": "BALSAMIC/config/reference.json",
        "sample": "BALSAMIC/config/sample.json",
        "analysis_paired": "BALSAMIC/config/analysis_paired.json",
        "analysis_paired_umi": "BALSAMIC/config/analysis_paired_umi.json",
        "analysis_single": "BALSAMIC/config/analysis_single.json",
        "analysis_single_umi": "BALSAMIC/config/analysis_single_umi.json",
        "panel_bed_file": "tests/test_data/references/GRCh37/panel/panel.bed",
        "test_reference": "tests/test_data/references/reference.json"
    }


@pytest.fixture(scope='session')
def conda_yaml():
    """
    conda env config file paths
    """
    return {
        "balsamic": "BALSAMIC/conda_yaml/BALSAMIC.yaml",
        "cancer-core": "BALSAMIC/conda_yaml/Cancer-Core.yaml",
        "cancer-p27": "BALSAMIC/conda_yaml/Cancer-p27.yaml",
        "cancer-py36": "BALSAMIC/conda_yaml/Cacner-py36.yaml",
        "cancer-vardict": "BALSAMIC/conda_yaml/Cancer-vardict.yaml",
        "cancer-vep": "BALSAMIC/conda_yaml/Cancer-vep.yaml",
        "cancer-vt": "BALSAMIC/conda_yaml/Cancer-vt.yaml"
    }
