import re
import subprocess
import json
from unittest import mock
from pathlib import Path


def test_run_analysis_dragen(invoke_cli, tumor_only_wgs_config):
    # GIVEN a WGS config file
    # WHEN running analysis
    result = invoke_cli(
        ['run', 'analysis', '-s', tumor_only_wgs_config, '--dragen'])

    # THEN it should run without any error
    assert result.exit_code == 0


def test_run_analysis_disable_variant_caller(invoke_cli, tumor_only_config):
    # GIVEN a tumor-only config file and variant caller to disable
    disabled_varcaller = "mutect"

    # WHEN running analysis
    result = invoke_cli([
        'run', 'analysis', '-s', tumor_only_config, '--disable-variant-caller',
        disabled_varcaller
    ])

    # THEN it should run without any error
    assert result.exit_code == 0
    assert disabled_varcaller not in result.output


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
    result = invoke_cli(
        ['run', 'analysis', '-s', tumor_only_config, '--mail-type', 'FAIL'])

    # THEN it should run without any error
    assert result.exit_code == 0


def test_run_analysis_click_abort(invoke_cli, tumor_only_config,
                                  tumor_normal_config):
    # GIVEN a tumor-only config file
    # WHEN running analysis
    result = invoke_cli(['run', 'analysis', '-r', '-s', tumor_only_config])

    # THEN it should abort with error
    assert result.exit_code == 1


def test_run_analysis_create_dir(invoke_cli, tumor_only_config):
    # GIVEN a tumor-only config file
    # WHEN running analysis

    with open(tumor_only_config) as fh:
        tumor_config = json.load(fh)
    log_dir = tumor_config['analysis']['log']
    Path(log_dir).mkdir(exist_ok=True)
    Path(log_dir, "logfile.log").touch(exist_ok=True)

    with mock.patch.object(subprocess, 'run') as mocked:
        mocked.return_value.stdout = 1
        invoke_cli([
            'run', 'analysis', '-s', tumor_only_config, '-r', '--account',
            'development'
        ])
        # THEN it should abort with error
        assert Path(re.sub('/$', '.1/', log_dir)).exists()
