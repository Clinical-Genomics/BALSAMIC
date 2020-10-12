import graphviz
from pathlib import Path
from unittest import mock
from BALSAMIC import __version__ as balsamic_version


def test_init_reference_write_json(invoke_cli, tmp_path,
                                   singularity_container):
    # Given test_reference.json
    test_genome_version = "hg19"
    test_new_dir = tmp_path / "test_reference_dir"
    test_new_dir.mkdir()

    # WHEN creating config.json in reference dir
    test_output_reference_config = test_new_dir / balsamic_version / test_genome_version / "config.json"
    test_output_reference_pdf = test_new_dir / balsamic_version / test_genome_version / "generate_ref_worflow_graph.pdf"

    result = invoke_cli([
        'init', 'reference', '-c', 'secret_key', '--singularity',
        singularity_container, '-o',
        str(test_new_dir)
    ])

    # THEN output config and pdf file generate and command exit code 0
    assert result.exit_code == 0
    assert Path(test_output_reference_pdf).exists()
    assert Path(test_output_reference_config).exists()


def test_init_reference_no_write_perm(
        tmp_path, invoke_cli, singularity_container, no_write_perm_path):
    # Given a path with no write permission
    test_new_dir = str(no_write_perm_path)

    # WHEN invoking config sample
    result = invoke_cli([
        'init', 'reference', '-c', 'secret_key', '--singularity',
        singularity_container, '-o',
        str(test_new_dir)
    ])

    # THEN it should create test_reference.json and exist with no error
    assert result.exit_code == 1


def test_init_reference_exception(invoke_cli, tmp_path, singularity_container):
    # Given test_reference.json
    test_new_dir = tmp_path / "test_reference_dir"
    test_new_dir.mkdir()

    with mock.patch.object(graphviz, 'Source') as mocked:
        mocked.return_value = None
        result = invoke_cli([
            'init', 'reference', '-c', 'secret_key', '--singularity',
            singularity_container, '-o',
            str(test_new_dir)
        ])

    assert result.exit_code == 1


def test_init_reference(invoke_cli):
    # WHEN invoking run reference command
    result = invoke_cli(['init', 'reference', '--help'])

    # THEN It should show the help message with all params
    assert "--snakefile" in result.output
    assert "--cosmic-key" in result.output
    assert "--singularity" in result.output
    assert result.exit_code == 0
