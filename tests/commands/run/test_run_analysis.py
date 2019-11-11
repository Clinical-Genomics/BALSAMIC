import os
import re
import subprocess
import json
import click
import pytest
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


def test_run_analysis_creat_dir(invoke_cli, tumor_only_config, tmp_path):
    # GIVEN a tumor-only config file
    # WHEN running analysis
    with open(tumor_only_config) as fh:
        tumor_config = json.load(fh)
    log_dir = tumor_config['analysis']['log']
    log_file = os.path.join(log_dir, 'log.file')
    result = invoke_cli(['run', 'analysis', '-s', tumor_only_config, '-r', '--account',
                         'development'])
    subprocess.run('touch ' + log_file, shell=True)

    with mock.patch.object(subprocess, 'run') as mocked:
        mocked.return_value.stdout = 1
        result = invoke_cli(['run', 'analysis', '-s', tumor_only_config, '-r', '--account',
                             'development'])
    # THEN it should abort with error
    assert os.path.exists(log_dir)
    assert os.path.exists(re.sub('/$', '.1/', log_dir))


def test_run_analysis_exception(invoke_cli, tumor_only_config):
    # GIVEN a tumor-only config file
    # WHEN running analysis with dummy option
    result = invoke_cli(['run', 'analysis', '-s', tumor_only_config, '-r', '--account',
                         'development', '--snakemake-opt', '--dummy'])
    # THEN It should abort the analysis with exit_code 1
    assert result.exit_code == 1
