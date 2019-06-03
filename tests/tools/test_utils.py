import unittest
import json

from BALSAMIC.tools import get_ref_path, iterdict
from BALSAMIC.tools.cli_utils import SnakeMake


class TestUtils(unittest.TestCase):
    """docstring for TestUtils"""

    def setUp(self):
        self.test_ref_json_path = "tests/test_data/references/reference.json"

    def test_get_ref_path(self):
        # GIVEN a sample json file path
        test_ref = self.test_ref_json_path

        # WHEN giving a path for json file,
        test_ref_json = get_ref_path(test_ref)

        # THEN It will read the file and return a dict with updated absolute path
        assert isinstance(test_ref_json, dict)
        assert test_ref_json['path']['genomefa'].startswith('/')

    def test_iterdict(self):
        """ GIVEN a dict for iteration """
        test_dict = json.load(open(self.test_ref_json_path, 'r'))

        # WHEN passing dict to this function
        dict_gen = iterdict(test_dict)

        # THEN it will create dict generator, we can iterate it, get the key, values as string
        for key, value in dict_gen:
            assert isinstance(key, str)
            assert isinstance(value, str)

    def test_snakemake_local(self):
        # GIVEN required params
        snakemake_local = SnakeMake()
        snakemake_local.working_dir = "/tmp/snakemake"
        snakemake_local.snakefile = "worflow/variantCalling_paired"
        snakemake_local.configfile = "sample_config.json"
        snakemake_local.run_mode = "local"

        # WHEN calling the build command
        shell_command = snakemake_local.build_cmd()

        # THEN it will contruct the snakemake command to run
        assert isinstance(shell_command, str)
        assert "worflow/variantCalling_paired" in shell_command
        assert "sample_config.json" in shell_command
        assert "/tmp/snakemake" in shell_command
        assert "--dryrun" in shell_command

    def test_snakemake_slurm(self):
        # GIVEN required params
        snakemake_slurm = SnakeMake()
        snakemake_slurm.sample_name = "test_sample"
        snakemake_slurm.working_dir = "/tmp/snakemake"
        snakemake_slurm.snakefile = "worflow/variantCalling_paired"
        snakemake_slurm.configfile = "sample_config.json"
        snakemake_slurm.run_mode = "slurm"
        snakemake_slurm.cluster_config = "cluster_config.json"
        snakemake_slurm.sbatch_py = "runsbatch.py"
        snakemake_slurm.log_path = "logs/"
        snakemake_slurm.script_path = "scripts/"
        snakemake_slurm.result_path = "results/"
        snakemake_slurm.qos = "normal"
        snakemake_slurm.sm_opt = ("containers")
        snakemake_slurm.run_analysis = True

        # WHEN calling the build command
        shell_command = snakemake_slurm.build_cmd()

        #THEN constructing snakecommand for slurm runner
        assert isinstance(shell_command, str)
        assert "worflow/variantCalling_paired" in shell_command
        assert "sample_config.json" in shell_command
        assert "/tmp/snakemake" in shell_command
        assert "--dryrun" not in shell_command
        assert "runsbatch.py" in shell_command
        assert "test_sample" in shell_command
        assert "containers" in shell_command
