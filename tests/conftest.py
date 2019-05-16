import pytest

from functools import partial
from click.testing import CliRunner

from BALSAMIC.commands import cli


@pytest.fixture
def cli_runner():
    runner = CliRunner()
    return runner


@pytest.fixture
def invoke_cli(cli_runner):
    return partial(cli_runner.invoke, cli)


@pytest.fixture(scope='session')
def config_files():
    return {
        "install": "BALSAMIC/config/install.json",
        "reference": "BALSAMIC/config/reference.json",
        "sample": "BALSAMIC/config/sample.json",
        "analysis_paired": "BALSAMIC/config/analysis_paired.json",
        "analysis_paired_umi": "BALSAMIC/config/analysis_paired_umi.json",
        "analysis_single": "BALSAMIC/config/analysis_single.json",
        "analysis_single_umi": "BALSAMIC/config/analysis_single_umi.json"
    }
