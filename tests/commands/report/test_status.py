import pytest

import BALSAMIC
from BALSAMIC.commands.base import cli


def test_status_tumor_only_panel(invoke_cli, tumor_only_config):
    # GIVEN a tumor-normal config file
    # WHEN running analysis
    result = invoke_cli(
        ['report', 'status', '--sample-config', tumor_only_config])

    # THEN it should run without any error
    assert result.exit_code == 0


def test_status_tumor_normal_panel(invoke_cli, tumor_normal_config):
    # GIVEN a tumor-normal config file
    # WHEN running analysis
    result = invoke_cli(
        ['report', 'status', '--sample-config', tumor_normal_config])

    # THEN it should run without any error
    assert result.exit_code == 0
