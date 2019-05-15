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
