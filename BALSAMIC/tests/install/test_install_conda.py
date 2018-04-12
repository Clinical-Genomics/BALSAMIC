from BALSAMIC.install import conda_env_check, conda_default_prefix, get_prefix
from unittest import TestCase, mock
#from unittest.mock import mock_open
import subprocess
import yaml


def test_conda_env_check():
    with mock.patch.object(subprocess, 'check_output') as mocked:
        mocked.return_value = '{"envs": ["/path/to/env/prefix_1", "/path/to/env/prefix_2" ]}'
        env_prefix = '/path/to/env/prefix_1'
        assert conda_env_check(env_prefix) == True


def test_conda_default_prefix():
    with mock.patch.object(subprocess, 'check_output') as mocked:
        mocked.return_value = '{"envs": ["/path/to/env/prefix_1", "/path/to/env/prefix_2" ], "conda_prefix": "default_prefix"}'
        env_prefix = 'default_prefix'
        assert conda_default_prefix() == env_prefix


def test_get_prefix():
    env_yaml = yaml.load(
        '{"channels": ["chan1", "chan2"], "dependencies": ["dep1", "dep2"], "prefix": "prefix1"}'
    )
    with mock.patch.object(yaml, 'load') as mocked_yaml:
        mocked_yaml.return_value = env_yaml
        assert get_prefix("mock_yaml_file") == "prefix1"
        assert not get_prefix("mock_yaml_file") == "chan1"

    env_yaml = yaml.load(
        '{"channels": ["chan1", "chan2"], "dependencies": ["dep1", "dep2"]}')
    with mock.patch.object(yaml, 'load') as mocked_yaml:
        mocked_yaml.return_value = env_yaml
        assert get_prefix("mock_yaml_file") == False
