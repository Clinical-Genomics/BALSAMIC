import BALSAMIC


def test_cli(invoke_cli):
    """
	BALSAMIC command cli testing
	
	by invoking --version 
	"""
    result = invoke_cli(['--version'])

    assert BALSAMIC.__version__ in result.output
