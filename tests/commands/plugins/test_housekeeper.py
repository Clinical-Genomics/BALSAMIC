import pytest

import BALSAMIC
from BALSAMIC.commands.base import cli

def test_housekeeper(invoke_cli, tumor_normal_config):
    # GIVEN a tumor-normal config file
    # WHEN running analysis
    result = invoke_cli(['plugins', 'housekeeper', '--sample-config', tumor_normal_config])

    # THEN it should run without any error
    assert result.exit_code == 0
