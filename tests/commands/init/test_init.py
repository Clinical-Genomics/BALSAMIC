import subprocess
import logging
import graphviz

from pathlib import Path
from unittest import mock
from BALSAMIC import __version__ as balsamic_version


def test_init_reference_write_json(invoke_cli, tmp_path,
                                   ):
    # Given test_reference.json
    test_genome_version = "hg19"
    test_container_version = "develop"
    test_new_dir = tmp_path / "test_reference_dir"
    test_new_dir.mkdir()

    # WHEN creating config.json in reference dir
    test_output_reference_config = test_new_dir / balsamic_version / test_genome_version / "config.json"
    test_output_reference_pdf = test_new_dir / balsamic_version / test_genome_version / "generate_ref_worflow_graph.pdf"

    result = invoke_cli([
        'init',
        '-o',
        str(test_new_dir),
        '-c',
        'secret_key',
        '-v',
        test_container_version,
    ])

    # THEN output config and pdf file generate and command exit code 0
    assert result.exit_code == 0
    assert Path(test_output_reference_pdf).exists()
    assert Path(test_output_reference_config).exists()


def test_init_reference_no_write_perm(
        tmp_path, invoke_cli, no_write_perm_path):
    # Given a path with no write permission
    test_genome_version = "hg19"
    test_container_version = "develop"
    test_new_dir = str(no_write_perm_path)

    # WHEN invoking config sample
    result = invoke_cli([
        'init',
        '-o',
        str(test_new_dir),
        '-c',
        'secret_key',
        '-v',
        test_container_version,
        '-g',
        test_genome_version,
    ])

    # THEN it should create test_reference.json and exist with no error
    assert result.exit_code == 1


def test_init_reference_graph_exception(invoke_cli, tmp_path):
    # Given test_reference.json
    test_new_dir = tmp_path / "test_reference_nonfunctional_graph"
    test_new_dir.mkdir()

    with mock.patch.object(graphviz, 'Source') as mocked:
        mocked.return_value = None
        result = invoke_cli([
            'init',
            '-o',
            str(test_new_dir),
            '-c',
            'secret_key',
        ])

    assert result.exit_code == 1


def test_init_container_force_dry(invoke_cli, tmp_path):
    # Given a dummy path
    test_new_dir = tmp_path / "test_container_dry_force"
    test_new_dir.mkdir()
    test_container_version = "develop"

    # WHEN force pull dry-run container
    result = invoke_cli([
        'init',
        '--outdir',
        str(test_new_dir),
        '-c',
        'secret_key',
        '--force',
        '-v',
        test_container_version,
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
            '-c',
            'secret_key',
        '--container-version',
        dummy_tag,
    ])

    # THEN command exit code 0
    assert result.exit_code == 0


def test_init_container_without_dry_run(invoke_cli, tmp_path):
    # Given a dummy path
    test_new_dir = tmp_path / "test_container_dir"
    test_new_dir.mkdir()

    with mock.patch.object(subprocess, 'run') as mocked:
        mocked.return_value = 0

        # WHEN pulling a container in a non dry-run mode
        result = invoke_cli([
            'init',
            '--outdir',
            str(test_new_dir),
            '-c',
            'secret_key',
            '--run-analysis',
        ])

        # THEN output config and pdf file generate and command exit code 0
        assert result.exit_code == 0


def test_init_container_wrong_tag(invoke_cli, tmp_path, caplog):
    # Given a dummy path
    test_new_dir = tmp_path / "test_container_dir"
    test_new_dir.mkdir()
    dummy_tag = "some_tag_that_does_not_exist_ngrtf123jsds3wqe2"

    # WHEN pulling a wrong container tag
    result = invoke_cli([
        'init',
        '--outdir',
        str(test_new_dir),
            '-c',
            'secret_key',
        '--container-version',
        dummy_tag,
    ])

    # THEN capture error log and error code
    assert result.exit_code > 0
