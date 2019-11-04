import pytest
import glob
from pathlib import Path
from click.testing import CliRunner

import BALSAMIC
from BALSAMIC.commands.base import cli


def test_run_analysis_tumor_normal_dry_run(invoke_cli, tumor_normal_config):
    # GIVEN a tumor-normal config file
    # WHEN running analysis
    result = invoke_cli(['run', 'analysis', '-s', tumor_normal_config])

    # THEN it should run without any error
    assert result.exit_code == 0


def test_run_analysis_tumor_only_dry_run(invoke_cli, tumor_only_config,
                                         tumor_normal_config):
    # GIVEN a tumor-only config file
    # WHEN running analysis
    result = invoke_cli(['run', 'analysis', '-s', tumor_only_config])

    # THEN it should run without any error
    assert result.exit_code == 0
