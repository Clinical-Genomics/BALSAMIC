import subprocess
import logging
from unittest import mock


def test_init_container(invoke_cli, tmp_path):
    # Given a dummy path
    test_new_dir = tmp_path / "test_container_dir"
    test_new_dir.mkdir()

    # WHEN when pulling container as dry-run
    result = invoke_cli([
        'init',
        '--outdir',
        str(test_new_dir),
        'container',
        '--dry',
    ])

    # THEN command exit code 0
    assert result.exit_code == 0


def test_init_container_force(invoke_cli, tmp_path):
    # Given a dummy path
    test_new_dir = tmp_path / "test_container_dir"
    test_new_dir.mkdir()

    # WHEN force pull dry-run container
    result = invoke_cli([
        'init',
        '--outdir',
        str(test_new_dir),
        'container',
        '--force',
        '--dry',
    ])

    # THEN command exit code 0
    assert result.exit_code == 0


def test_init_container_specific_tag(invoke_cli, tmp_path):
    # Given a dummy path
    test_new_dir = tmp_path / "test_container_dir"
    test_new_dir.mkdir()
    dummy_tag = "develop"

    # WHEN pulling a specific tag other than standard version
    result = invoke_cli([
        'init',
        '--outdir',
        str(test_new_dir),
        'container',
        '--force',
        '--dry',
        '--container-version',
        dummy_tag,
    ])

    # THEN command exit code 0
    assert result.exit_code == 0


def test_init_container_without_dry_run(invoke_cli, tmp_path):
    # Given a dummy path
    test_new_dir = tmp_path / "test_container_dir"
    test_new_dir.mkdir()

    with mock.patch.object(subprocess, 'check_output') as mocked:
        mocked.return_value = 0

        # WHEN pulling a container in a non dry-run mode
        result = invoke_cli([
            'init',
            '--outdir',
            str(test_new_dir),
            'container',
        ])

        # THEN output config and pdf file generate and command exit code 0
        assert result.exit_code == 0


def test_init_container_capture_failed_download(invoke_cli, tmp_path, caplog):
    # Given a dummy path
    test_new_dir = tmp_path / "test_container_dir"
    test_new_dir.mkdir()
    dummy_tag = "develop"

    with mock.patch.object(subprocess,
                           'check_output') as mocked, caplog.at_level(
                               logging.ERROR):
        mocked.side_effect = 1

        # WHEN pulling a wrong container tag
        result = invoke_cli([
            'init',
            '--outdir',
            str(test_new_dir),
            'container',
            '--container-version',
            dummy_tag,
        ])

        # THEN capture error log and error code
        assert "Failed to pull singularity image" in caplog.text
        assert result.exit_code > 0
