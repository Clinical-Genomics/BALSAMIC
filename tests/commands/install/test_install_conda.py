from BALSAMIC.commands.install import conda_env_check, conda_default_prefix, get_prefix
from unittest import TestCase, mock
import subprocess
import yaml


def test_conda_env_check():
    with mock.patch.object(subprocess, 'check_output') as mocked:
        # GIVEN a dict of environment prefixes installed
        mocked.return_value = '{"envs": ["/path/to/env/prefix_1", "/path/to/env/prefix_2" ]}'

        # WHEN matching an existing conda environment
        env_prefix = '/path/to/env/prefix_1'
        actual = conda_env_check(env_prefix)

        # THEN assert that environment exists
        assert actual == True


def test_conda_default_prefix():
    with mock.patch.object(subprocess, 'check_output') as mocked:
        # GIVEN a dict of conda environment config output
        mocked.return_value = '{"envs": ["/path/to/env/prefix_1", "/path/to/env/prefix_2" ], "conda_prefix": "prefix_1"}'

        # WHEN get default environment prefix
        actual = conda_default_prefix()

        # THEN assert if the default_prefix matches correctly
        env_prefix = 'prefix_1'
        assert actual == env_prefix


def test_get_prefix_with_prefix_key():

    # GIVEN conda env yaml with prefix entry
    env_yaml = yaml.load(
        '{"channels": ["chan1", "chan2"], "dependencies": ["dep1", "dep2"], "prefix": "prefix1"}'
    )

    # WHEN yaml read and returned
    with mock.patch.object(yaml, 'load') as mocked_yaml:
        mocked_yaml.return_value = env_yaml

        # THEN assert if the yaml has "prefix1" value for prefix entry
        assert get_prefix("mock_yaml_file") == "prefix1"
        assert not get_prefix("mock_yaml_file") == "chan1"

def test_get_prefix_without_prefix_key():
    # GIVEN conda env yaml without prefix entry
    env_yaml = yaml.load(
        '{"channels": ["chan1", "chan2"], "dependencies": ["dep1", "dep2"]}')

    # WHEN yaml read and return
    with mock.patch.object(yaml, 'load') as mocked_yaml:
        mocked_yaml.return_value = env_yaml

        # THEN make sure the output is False
        assert get_prefix("mock_yaml_file") == False
