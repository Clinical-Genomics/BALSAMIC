from BALSAMIC.install import conda_env_check
import json
import subprocess

#def test_conda_env_check(conda_env_json, _env_prefix)
#  assert existing_conda_env == True


def test_conda_env_check():
    conda_env_json = json.loads(subprocess.check_output(["conda", "env", "list", "--json"], stderr=subprocess.STDOUT))
    
    p = conda_env_json["envs"][0]
    assert conda_env_check(p) == True 


#def conda_env_check(env_prefix):
#    """
#    Check if conda env exists.
#    Input: Conda env prefix extracted from conda yaml file
#    Output: True/False
#    """
#    try:
#        p = json.loads(
#            subprocess.check_output(
#                ["conda", "env", "list", "--json"], stderr=subprocess.STDOUT))
#    except subprocess.CalledProcessError as e:
#        print(e.output.decode())
#        if verbose:
#            raise e.output.decode()
#
#    return env_prefix in p["envs"]
