import os
import re
import subprocess
import json
import click
import pytest
import glob
from unittest import mock
from pathlib import Path


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


def test_run_analysis_create_dir(tmpdir_factory, invoke_cli, tumor_only_config, analysis_dir, tmp_path):
    # GIVEN a tumor-only config file
    # WHEN running analysis
    case_name = Path(tumor_only_config).name.split(".")[0]
    dummy_log = tmpdir_factory.mktemp(
        analysis_dir + "/" + case_name + '/logs',
        numbered=False).join('example_file').write('not_empty')
    dummy_balsamic_stat = tmpdir_factory.mktemp(
        analysis_dir + "/" + case_name + '/analysis/vep',
        numbered=False).join('vcf_merge.balsamic_stat').write('not_empty')

    with open(tumor_only_config) as fh:
        tumor_config = json.load(fh)
    log_dir = tumor_config['analysis']['log']

    with mock.patch.object(subprocess, 'run') as mocked:
        mocked.return_value.stdout = 1
        invoke_cli(['run', 'analysis', '-s', tumor_only_config, '-r', '--account',
                             'development'])
        # THEN it should abort with error
        assert Path(re.sub('/$', '.1/', log_dir)).exists()


# def test_run_analysis_exception(invoke_cli, tumor_only_config):
#     # GIVEN a tumor-only config file
#     # WHEN running analysis with dummy option
#     with mock.patch.object(subprocess, 'run') as mocked:
#         result = invoke_cli(['run', 'analysis', '-s', tumor_only_config, '-r', '--account',
#                              'development', '--snakemake-opt', '--dummy'])
#     # THEN It should abort the analysis with exit_code 1
#     assert result.exit_code == 1
