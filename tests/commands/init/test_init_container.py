import subprocess
import logging
from unittest import mock

import click
import pytest

from BALSAMIC.utils.exc import BalsamicError


def test_init_container(invoke_cli, tmp_path):
    # Given a dummy path
    test_new_dir = tmp_path / "test_container_dir"
    test_new_dir.mkdir()

    # WHEN creating config.json in reference dir
    result = invoke_cli(
        ['init', 'container', '--dry', '--out-dir',
         str(test_new_dir)])

    # THEN output config and pdf file generate and command exit code 0
    assert result.exit_code == 0


def test_init_container_force(invoke_cli, tmp_path):
    # Given a dummy path
    test_new_dir = tmp_path / "test_container_dir"
    test_new_dir.mkdir()

    # WHEN creating config.json in reference dir
    result = invoke_cli([
        'init', 'container', '--force', '--dry', '--out-dir',
        str(test_new_dir)
    ])

    # THEN output config and pdf file generate and command exit code 0
    assert result.exit_code == 0


def test_init_container_specific_tag(invoke_cli, tmp_path):
    # Given a dummy path
    test_new_dir = tmp_path / "test_container_dir"
    test_new_dir.mkdir()
    dummy_tag = "cool_new_feature"

    # WHEN creating config.json in reference dir
    result = invoke_cli([
        'init', 'container', '--force', '--dry', '--container-version',
        dummy_tag, '--out-dir',
        str(test_new_dir)
    ])

    # THEN output config and pdf file generate and command exit code 0
    assert result.exit_code == 0


def test_init_container_without_dry_run(invoke_cli, tmp_path):
    # Given a dummy path
    test_new_dir = tmp_path / "test_container_dir"
    test_new_dir.mkdir()
    dummy_tag = "cool_new_feature"

    # WHEN calling scheduler_main with mocked subprocess
    with mock.patch.object(subprocess, 'check_output') as mocked:
        # WHEN creating config.json in reference dir
        result = invoke_cli([
            'init', 'container', '--force', '--container-version', dummy_tag,
            '--out-dir',
            str(test_new_dir)
        ])

        # THEN output config and pdf file generate and command exit code 0
        assert result.exit_code == 0 


def test_init_container_capture_failed_download(invoke_cli, tmp_path, caplog):
    # Given a dummy path
    test_new_dir = tmp_path / "test_container_dir"
    test_new_dir.mkdir()
    dummy_tag = "some_tag_that_does_not_exist_ngrtf123jsds3wqe2"

    # WHEN calling scheduler_main with mocked subprocess
    with caplog.at_level(logging.ERROR):
        # WHEN creating config.json in reference dir
        result = invoke_cli([
            'init', 'container', '--container-version', dummy_tag, 
            '--out-dir',
            str(test_new_dir)
        ])
      
        assert "Failed to pull singularity image" in caplog.text 
