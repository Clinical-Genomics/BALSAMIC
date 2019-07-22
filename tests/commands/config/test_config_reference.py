import pytest
import glob
from pathlib import Path
from click.testing import CliRunner

from BALSAMIC.commands.base import cli


def test_config_reference_write_json(install_config, invoke_cli, tmp_path):
    # Given test_reference.json
    test_new_dir = tmp_path / "test_reference_dir"
    test_new_dir.mkdir()

    # WHEN creating config.json in reference dir
    test_output_reference_config = test_new_dir / "config.json"
    test_output_reference_pdf = test_new_dir / "generate_ref_worflow_graph.pdf"

    result = invoke_cli([
        'config', 'reference', '-i', install_config, '-c', 'secret_key', '-o',
        str(test_new_dir)
    ])

    # THEN output config and pdf file generate and command exit code 0
    assert result.exit_code == 0
    assert Path(test_output_reference_pdf).exists()
    assert Path(test_output_reference_config).exists()


def test_config_reference_no_write_perm(tmp_path, invoke_cli,
                                        no_write_perm_path):
    # Given a path with no write permission
    test_new_dir = str(no_write_perm_path)

    # WHEN invoking config sample
    result = invoke_cli(
        ['config', 'reference', '-c', 'secret_key', '-o',
         str(test_new_dir)])

    # THEN it should create test_reference.json and exist with no error
    assert result.exit_code == 1

