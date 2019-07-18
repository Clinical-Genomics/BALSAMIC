from pathlib import Path
import BALSAMIC
import pytest


def test_cli(invoke_cli):
    # GIVEN I want to see version of the program
    # WHEN I am asking to see version
    result = invoke_cli(['--version'])

    # THEN It should show the version of the program
    assert BALSAMIC.__version__ in result.output


def test_config(invoke_cli):
    # GIVEN I want to see config command options
    # WHEN asking to show config options
    result = invoke_cli(['config'])

    # THEN It should show config options in result
    assert 'sample' in result.output


def test_sample(invoke_cli):
    # GIVEN want to see config-sample params with help option
    # WHEN asking to show params for config-sample
    result = invoke_cli(['config', 'sample', '--help'])

    # THEN It should show all params reuired for config-sample
    assert 'sample-id' in result.output
    assert result.exit_code == 0


def test_sample_missing_opt(invoke_cli):
    # WHEN invoking command with missing options
    result = invoke_cli(['config', 'sample'])

    # THEN It should throw missiong option error
    assert 'Error: Missing option' in result.output
    assert result.exit_code == 2


def test_install(invoke_cli):
    # WHEN invoking install command with help option
    result = invoke_cli(['install', '--help'])

    # THEN It should list all the params without error
    assert result.exit_code == 0
    assert '--input-conda-yaml' in result.output


def test_install_missing_opt(invoke_cli):
    # WHEN invoking install command with missing option
    result = invoke_cli(['install'])

    # THEN It should esclate the missing option error
    assert 'Error: Missing option' in result.output
    assert result.exit_code == 2


def test_install_invalid(invoke_cli):
    # GIVEN wrong wrong env type, '-t' should have 'P','D','S'.
    # WHEN invoking command with invalid input
    result = invoke_cli(['install', '-t', 'foo'])

    # THEN It should throw invalid input error
    assert result.exit_code == 2
    assert 'Error: Invalid value' in result.output


def test_run(invoke_cli):
    # WHEN asking to options for run command
    result = invoke_cli(['run', '--help'])

    # THEN It should show all the params without any error
    assert result.exit_code == 0
    assert "analysis" in result.output


def test_run_analysis(invoke_cli):
    # WHEN invoking run analysis command
    result = invoke_cli(['run', 'analysis', '--help'])

    # THEN it should show all params without error
    assert "--analysis-type" in result.output
    assert "--snake-file" in result.output
    assert "--sample-config" in result.output
    assert "--run-mode" in result.output
    assert "--cluster-config" in result.output
    assert "--run-analysis" in result.output


def test_run_missing_opt(invoke_cli):
    # WHEN invoking run command with missing option
    result = invoke_cli(['run', 'analysis'])

    # THEN It should throw missiong option error
    assert 'Error: Missing option' in result.output
    assert result.exit_code == 2


def test_run_analysis_invalid(invoke_cli):
    # WHEN invoking run with invalid input value
    result = invoke_cli(['run', 'analysis', '--run-mode', 'foo'])

    # THEN It should throw invalid value error
    assert result.exit_code == 2
    assert 'Error: Invalid value' in result.output


def test_run_reference(invoke_cli):
    # WHEN invoking run reference command
    result = invoke_cli(['run', 'reference', '--help'])

    # THEN It should show the help message with all params
    assert "--snakefile" in result.output
    assert "--configfile" in result.output
    assert "--run-mode" in result.output
    assert "--cluster-config" in result.output
    assert "--run-analysis" in result.output
    assert result.exit_code == 0


def test_run_ref_invalid(invoke_cli):
    # WHEN invoking run reference command with invalid param
    result = invoke_cli(['run', 'reference', '--run-mode', 'foo'])

    # THEN It should throw invalid value error
    assert result.exit_code == 2
    assert 'Error: Invalid value' in result.output

def test_config_reference(tmp_path, invoke_cli):
    # Given test_reference.json
    test_new_dir = tmp_path / "test_reference_dir" 
    test_new_dir.mkdir()
    test_new_dir.chmod(0o777)

    test_output_reference_config = test_new_dir / "config.json" 
    test_output_reference_pdf = test_new_dir / "generate_ref_dag.pdf" 

    #test_output_reference_config.touch()
    #test_output_reference_pdf.touch()

    # WHEN invoking config sample
    result = invoke_cli(['config', 'reference', '-c', 'secret_key', '-o', str(test_new_dir)])

    # THEN it should create test_reference.json and exist with no error
    #assert result.exit_code == 0
    import glob
    print(glob.glog(str(test_new_dir)+"/*"))
    assert Path(str(test_output_reference_config)).exists()
    assert Path(str(test_output_reference_pdf)).exists()

def test_config_reference_no_write_perm(tmp_path, invoke_cli, no_write_perm_path):
    # Given a path with no write permission
    test_new_dir = str(no_write_perm_path)

    # WHEN invoking config sample
    result = invoke_cli(['config', 'reference', '-c', 'secret_key', '-o', str(test_new_dir)])

    # THEN it should create test_reference.json and exist with no error
    assert result.exit_code == 1
