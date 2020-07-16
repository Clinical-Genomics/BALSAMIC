import pytest

import BALSAMIC
from BALSAMIC.commands.base import cli


def test_scout_tumor_normal(invoke_cli, tumor_normal_config):
    # GIVEN a tumor-normal config file
    # WHEN running analysis
    result = invoke_cli([
        'plugins', 'scout', '--sample-config', tumor_normal_config,
        '--customer-id', 'cust000'
    ])

    # THEN it should run without any error
    print(result)
    assert result.exit_code == 0


def test_scout_tumor_only(invoke_cli, tumor_only_config):
    # GIVEN a tumor-only config file
    # WHEN running analysis
    result = invoke_cli([
        'plugins', 'scout', '--sample-config', tumor_only_config,
        '--customer-id', 'cust000'
    ])

    # THEN it should run without any error
    print(result)
    assert result.exit_code == 0
