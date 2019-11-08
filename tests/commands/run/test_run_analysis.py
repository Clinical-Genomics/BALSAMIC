import os
import subprocess
import json
from unittest import mock


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
    result = invoke_cli(['run', 'analysis', '-s', tumor_only_config, '--mail-type', 'FAIL'])

    # THEN it should run without any error
    assert result.exit_code == 0


def test_run_analysis_click_abort(invoke_cli, tumor_only_config, tumor_normal_config):
    # GIVEN a tumor-only config file
    # WHEN running analysis
    result = invoke_cli(['run', 'analysis', '-r', '-s', tumor_only_config])

    # THEN it should abort with error
    assert result.exit_code == 1


def test_run_analysis_creat_dir(invoke_cli, tumor_only_config):
    # GIVEN a tumor-only config file
    # WHEN running analysis
    with open(tumor_only_config) as fh:
        tumor_config = json.load(fh)
    scripts_dir = tumor_config['analysis']['script']
    with mock.patch.object(subprocess, 'run') as mocked:
        mocked.return_value.stdout = 1
        result = invoke_cli(['run', 'analysis', '-s', tumor_only_config, '-r'])
        result = invoke_cli(['run', 'analysis', '-s', tumor_only_config, '-r'])
    # THEN it should abort with error
    assert os.path.exists(scripts_dir)
    assert os.path.exists(scripts_dir.replace('$/', '.1/'))
