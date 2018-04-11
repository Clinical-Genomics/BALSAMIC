from BALSAMIC.install import conda_env_check, conda_default_prefix
from unittest import TestCase, mock
import subprocess
import json

def test_conda_env_check():
    with mock.patch.object(subprocess, 'check_output') as mocked:
        mocked.return_value = '{"envs": ["/path/to/env/prefix_1", "/path/to/env/prefix_2" ]}'
        env_prefix = "/path/to/env/prefix_1"
        assert conda_env_check(env_prefix) == True

def test_conda_default_prefix():
    with mock.patch.object(subprocess, 'check_output') as mocked:
        mocked.return_value = '{"envs": ["/path/to/env/prefix_1", "/path/to/env/prefix_2" ], "conda_prefix": "default_prefix"}'
        env_prefix = "default_prefix"
        assert conda_default_prefix() == env_prefix

#def conda_default_prefix():
#    try:
#        p = json.loads(
#            subprocess.check_output(
#                ["conda", "info", "--json"], stderr=subprocess.STDOUT))
#    except subprocess.CalledProcessError as e:
#        print(e.output.decode())
#        if verbose:
#            raise e.output.decode()
#
#    return p["conda_prefix"]
