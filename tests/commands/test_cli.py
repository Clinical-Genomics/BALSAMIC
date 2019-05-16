import BALSAMIC


def test_cli(invoke_cli):
    """
	BALSAMIC command cli testing
	
	by invoking --version 
	"""
    result = invoke_cli(['--version'])

    assert BALSAMIC.__version__ in result.output


def test_config(invoke_cli):
    """ config command testing - It returns two commands in usage"""
    result = invoke_cli(['config'])

    assert 'sample' in result.output
    assert 'report' in result.output


def test_sample(invoke_cli):
    """ balsamic config sample --help - returns usage with all required params """
    result = invoke_cli(['config', 'sample', '--help'])

    assert 'sample-id' in result.output
    assert result.exit_code == 0


def test_sample_missing_opt(invoke_cli):
    """ testing the missing option"""
    result = invoke_cli(['config', 'sample'])

    assert 'Error: Missing option' in result.output
    assert result.exit_code == 2


def test_sample_invalid(invoke_cli):
    """ Invalid option for analysis type """
    result = invoke_cli(['config', 'sample', '-a', 'foo'])

    assert result.exit_code == 2
    assert 'Error: Invalid' in result.output


def test_install(invoke_cli):

    result = invoke_cli(['install', '--help'])

    assert result.exit_code == 0
    assert '--input-conda-yaml' in result.output


def test_install_missing_opt(invoke_cli):

    result = invoke_cli(['install'])

    assert 'Error: Missing option' in result.output
    assert result.exit_code == 2


def test_install_invalid(invoke_cli):
    """ Invalid option for environment type """
    result = invoke_cli(['install', '-t', 'foo'])

    assert result.exit_code == 2
    assert 'Error: Invalid' in result.output
