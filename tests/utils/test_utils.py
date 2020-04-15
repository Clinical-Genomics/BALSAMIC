import json
import os
import pytest
import sys
import copy
import collections
from pathlib import Path

from BALSAMIC.utils.cli import SnakeMake
from BALSAMIC.utils.cli import CaptureStdout
from BALSAMIC.utils.cli import iterdict
from BALSAMIC.utils.cli import get_packages
from BALSAMIC.utils.cli import get_package_split
from BALSAMIC.utils.cli import get_snakefile
from BALSAMIC.utils.cli import createDir
from BALSAMIC.utils.cli import get_ref_path
from BALSAMIC.utils.cli import write_json
from BALSAMIC.utils.cli import get_config
from BALSAMIC.utils.cli import recursive_default_dict
from BALSAMIC.utils.cli import convert_defaultdict_to_regular_dict
from BALSAMIC.utils.cli import get_file_status_string
from BALSAMIC.utils.cli import get_from_two_key
from BALSAMIC.utils.cli import find_file_index
from BALSAMIC.utils.rule import get_chrom
from BALSAMIC.utils.rule import get_vcf
from BALSAMIC.utils.rule import get_sample_type
from BALSAMIC.utils.rule import get_conda_env
from BALSAMIC.utils.rule import get_picard_mrkdup
from BALSAMIC.utils.rule import get_script_path
from BALSAMIC.utils.rule import get_result_dir
from BALSAMIC.utils.rule import get_threads


def test_recursive_default_dict():
    # GIVEN a dictionary
    test_dict = recursive_default_dict()
    test_dict['key_1']['key_2'] = 'value_1'

    # WHEN it is recursively creates a default dictionary
    # THEN the output should be a dicitionary
    assert isinstance(test_dict, collections.defaultdict)
    assert 'key_2' in test_dict['key_1']


def test_convert_defaultdict_to_regular_dict():
    # GIVEN a recursively created default dict
    test_dict = recursive_default_dict()
    test_dict['key_1']['key_2'] = 'value_1'

    # WHEN converting it back to normal dict
    test_dict = convert_defaultdict_to_regular_dict(test_dict)
    
    # THEN the output type should be dict and not defaultdict
    assert not isinstance(test_dict, collections.defaultdict)
    assert isinstance(test_dict, dict)
    assert 'key_2' in test_dict['key_1']


def test_get_ref_path(config_files):
    # GIVEN a sample json file path
    test_ref = config_files['test_reference']

    # WHEN giving a path for json file,
    test_ref_json = get_ref_path(test_ref)

    # THEN It will read the file and return a dict with updated absolute path
    assert isinstance(test_ref_json, dict)
    for ref, ref_path in test_ref_json['reference'].items():
        assert Path(ref_path).exists()


def test_iterdict(config_files):
    """ GIVEN a dict for iteration """
    test_dict = json.load(open(config_files['test_reference'], 'r'))

    # WHEN passing dict to this function
    dict_gen = iterdict(test_dict)

    # THEN it will create dict generator, we can iterate it, get the key, values as string
    for key, value in dict_gen:
        assert isinstance(key, str)
        assert isinstance(value, str)


def test_snakemake_local():
    # GIVEN required params
    snakemake_local = SnakeMake()
    snakemake_local.working_dir = "/tmp/snakemake"
    snakemake_local.snakefile = "worflow/variantCalling_paired"
    snakemake_local.configfile = "sample_config.json"
    snakemake_local.run_mode = "local"
    snakemake_local.use_singularity = True
    snakemake_local.singularity_bind = ["path_1", "path_2"]
    snakemake_local.forceall = True

    # WHEN calling the build command
    shell_command = snakemake_local.build_cmd()

    # THEN it will contruct the snakemake command to run
    assert isinstance(shell_command, str)
    assert "worflow/variantCalling_paired" in shell_command
    assert "sample_config.json" in shell_command
    assert "/tmp/snakemake" in shell_command
    assert "--dryrun" in shell_command
    assert "--forceall" in shell_command


def test_snakemake_slurm():
    # GIVEN required params
    snakemake_slurm = SnakeMake()
    snakemake_slurm.case_name = "test_case"
    snakemake_slurm.working_dir = "/tmp/snakemake"
    snakemake_slurm.snakefile = "worflow/variantCalling_paired"
    snakemake_slurm.configfile = "sample_config.json"
    snakemake_slurm.run_mode = "cluster"
    snakemake_slurm.cluster_config = "cluster_config.json"
    snakemake_slurm.scheduler = "sbatch.py"
    snakemake_slurm.log_path = "logs/"
    snakemake_slurm.script_path = "scripts/"
    snakemake_slurm.result_path = "results/"
    snakemake_slurm.qos = "normal"
    snakemake_slurm.account = "development"
    snakemake_slurm.profile = "slurm"
    snakemake_slurm.mail_type = "FAIL"
    snakemake_slurm.mail_user = "john.doe@example.com"
    snakemake_slurm.sm_opt = ("containers", )
    snakemake_slurm.use_singularity = True
    snakemake_slurm.singularity_bind = ["path_1", "path_2"]
    snakemake_slurm.run_analysis = True

    # WHEN calling the build command
    shell_command = snakemake_slurm.build_cmd()

    # print(shell_command)
    # THEN constructing snakecommand for slurm runner
    assert isinstance(shell_command, str)
    assert "worflow/variantCalling_paired" in shell_command
    assert "sample_config.json" in shell_command
    assert "/tmp/snakemake" in shell_command
    assert "--dryrun" not in shell_command
    assert "sbatch.py" in shell_command
    assert "test_case" in shell_command
    assert "containers" in shell_command


def test_get_packages(conda):
    # GIVEN a conda yaml file
    balsamic_yaml = conda['balsamic-base']

    # WHEN passing conda yaml file to get_packages
    packages = get_packages(balsamic_yaml)

    # THEN It should return all tools(packages) in that yaml file
    assert any("pip=9" in tool for tool in packages)
    assert any("python=3" in tool for tool in packages)


def test_get_pacakges_error(conda):
    with pytest.raises(Exception, match=r"No such file or directory"):
        # GIVEN a wrong path for config yaml file
        invalid_yaml = "/tmp/balsamic.yaml"

        # WHEN passing this invalid path into get_packages
        # THEN It shoud raise the file not found error
        assert get_packages(invalid_yaml)


def test_get_package_split(conda):
    # GIVEN a list of conda config yaml file paths
    yamls = conda.values()

    # WHEN giving this list of yamls to get_packaged split
    packages = get_package_split(yamls)

    # THEN It should return a dictionary with tools information
    assert isinstance(packages, dict)
    assert "bwa" in packages
    assert "bcftools" in packages
    assert "fastqc" in packages


def test_get_script_path():
    # GIVEN list of scripts
    custom_scripts = ["refseq_sql.awk"]

    # WHEN asking to get path for the script
    for script in custom_scripts:
        # THEN assert full path for script
        assert get_script_path(script).startswith('/')

        # THEN assert if script file exists
        assert Path(get_script_path(script)).is_file()


def test_get_snakefile():
    # GIVEN analysis_type for snakemake workflow
    workflow = [("paired", "wgs"), ("paired", "targeted"), ("single", "wgs"),
                ("single", "targeted"), ("qc", ""), ("generate_ref", "")]

    # WHEN asking to see snakefile for paired
    for analysis_type, sequencing_type in workflow:
        snakefile = get_snakefile(analysis_type, sequencing_type)
        pipeline = ''

        if sequencing_type == 'targeted':
            pipline = "BALSAMIC/workflows/VariantCalling"
        elif sequencing_type == 'wgs':
            pipeline = "BALSAMIC/workflows/VariantCalling_sentieon"
        elif analysis_type == 'qc':
            pipeline = "BALSAMIC/workflows/Alignment"
        elif analysis_type == 'generate_ref':
            pipeline = "BALSAMIC/workflows/GenerateRef"

        # THEN it should return the snakefile path
        # THEN assert file exists
        assert snakefile.startswith('/')
        assert pipeline in snakefile
        assert Path(snakefile).is_file()


def test_get_chrom(config_files):
    # Given a panel bed file
    bed_file = config_files["panel_bed_file"]
    actual_chrom = ['10','11','16','17','18','19','2','3','4','6','7','9','X']

    # WHEN passing this bed file
    test_chrom = get_chrom(bed_file)

    # THEN It should return list of chrom presents in that bed file
    assert set(actual_chrom) == set(test_chrom)


def test_get_vcf(sample_config):
    # GIVEN a sample_config dict, varinat callers list
    variant_callers = ['mutect', 'vardict', 'manta']

    # WHEN passing args to that function
    vcf_list = get_vcf(sample_config, variant_callers,
                       [sample_config["analysis"]["case_id"]])

    # THEN It should return the list of vcf file names
    assert any("mutect" in vcf_name for vcf_name in vcf_list)
    assert any("vardict" in vcf_name for vcf_name in vcf_list)
    assert any("manta" in vcf_name for vcf_name in vcf_list)


def test_get_sample_type(sample_config):
    # GIVEN a sample_config dict, bio_type as tumor
    bio_type = 'tumor'

    # WHEN calling get_sample_type with bio_type
    sample_id = get_sample_type(sample_config["samples"], bio_type)

    # THEN It should return the tumor samples id
    assert sample_id == ['S1_R']


def test_get_picard_mrkdup(sample_config):
    # WHEN passing sample_config
    picard_str = get_picard_mrkdup(sample_config)

    # THEN It will return the picard str as rmdup
    assert "mrkdup" == picard_str


def test_get_picard_mrkdup_rmdup(sample_config):
    # WHEN passing sample_config
    sample_config_rmdup = copy.deepcopy(sample_config)
    sample_config_rmdup["QC"]["picard_rmdup"] = True

    picard_str = get_picard_mrkdup(sample_config_rmdup)

    # THEN It will return the picard str as rmdup
    assert "rmdup" == picard_str

def test_createDir(tmp_path):
    # GIVEN a directory path
    # WHEN directory path is not yet created
    test_new_dir = tmp_path / "new_dir"

    # THEN it should create and return dir name
    test_new_dir_created = createDir(str(test_new_dir))
    assert test_new_dir_created == str(tmp_path / "new_dir")
    assert Path(test_new_dir_created).is_dir()

    # GIVEN a directory path
    test_log_dir = tmp_path / "existing_dir"

    # WHEN directory path exists
    test_log_dir.mkdir()

    # THEN it should return log_dir name incremented
    test_log_dir_created = createDir(str(test_log_dir), [])
    assert test_log_dir_created == str(tmp_path / "existing_dir.1")
    assert Path(test_log_dir_created).is_dir()

    # GIVEN a directory path
    test_log_dir = tmp_path / "existing_dir_with_dot.1"

    # WHEN directory path exists
    test_log_dir.mkdir()

    # THEN it should return log_dir name incremented
    test_log_dir_created = createDir(str(test_log_dir), [])
    assert test_log_dir_created == str(tmp_path / "existing_dir_with_dot.2")
    assert Path(test_log_dir_created).is_dir()


def test_get_result_dir(sample_config):
    # WHEN a sample_config dict
    # GIVEN a sample_config dict
    # THEN get_result_dir should return result directory
    assert get_result_dir(sample_config) == sample_config["analysis"]["result"]


def test_get_conda_env_found(tmp_path):
    # GIVEN a balsamic_env yaml
    balsamic_env = "BALSAMIC/config/balsamic_env.yaml"

    # WHEN passing pkg name with this yaml file
    conda_env = get_conda_env(balsamic_env, 'cnvkit')

    # THEN It should return the conda env which has that pkg
    assert conda_env == "BALSAMIC_py36"


def test_get_conda_env_not_found(tmp_path):
    # GIVEN a balsamic_env yaml
    balsamic_env = "BALSAMIC/config/balsamic_env.yaml"

    # WHEN passing pkg name with this yaml file
    # THEN It should return the conda env which has that pkg
    with pytest.raises(KeyError):
        get_conda_env(balsamic_env, 'unknown_package')


def test_capturestdout():
    # GIVEN a catpurestdout context
    test_stdout_message = 'Message to stdout'
    with CaptureStdout() as captured_stdout_message:
        print(test_stdout_message, file=sys.stdout)

    assert "".join(captured_stdout_message) == test_stdout_message


def test_get_config():
    # GIVEN the config files name
    config_files = ["sample", "analysis"]
    # WHEN passing file names
    for config_file in config_files:
        # THEN return the config files path
        assert get_config(config_file)


def test_get_config_wrong_config():
    # GIVEN the config files name
    config_file = 'non_existing_config'

    # WHEN passing file names
    # THEN return the config files path
    with pytest.raises(FileNotFoundError):
        assert get_config(config_file)


def test_write_json(tmp_path, config_files):
    # GIVEN a dict from sample json file (reference.json)
    ref_json = json.load(open(config_files['reference'], 'r'))

    tmp = tmp_path / "tmp"
    tmp.mkdir()
    output_json = tmp / "output.json"

    # WHEN passing dict and file name
    write_json(ref_json, output_json)
    output = output_json.read_text()

    # THEN It will create a json file with given dict
    for key, value in iterdict(ref_json):
        assert key in output
        assert value in output

    assert len(list(tmp.iterdir())) == 1


def test_write_json_error(tmp_path, config_files):
    with pytest.raises(Exception, match=r"Is a directory"):
        # GIVEN a invalid dict
        ref_json = {"path": "/tmp", "reference": ""}
        tmp = tmp_path / "tmp"
        tmp.mkdir()
        output_json = tmp / "/"

        # WHEN passing a invalid dict
        # THEN It will raise the error
        assert write_json(ref_json, output_json)


def test_get_threads(config_files):
    # GIVEN cluster config file and rule name
    cluster_config = json.load(open(config_files['cluster_json'], 'r'))
    rule_name = 'sentieon_align_sort'

    # WHEN passing cluster_config and rule_name
    # THEN It should return threads value '12'
    assert get_threads(cluster_config, rule_name)


def test_get_file_status_string_file_exists(tmpdir):
    # GIVEN an existing file and condition_str False
    file_exist = tmpdir.mkdir("temporary_path").join("file_exists")
    file_exist.write("dummy_file_content")
    
    # WHEN checking for file string
    result = get_file_status_string(str(file_exist))

    # THEN it should not return empty str
    assert "Found" in result[0].value_no_colors


def test_get_file_status_string_file_not_exist():
    # GIVEN an existing file and condition_str False
    file_not_exist = "some_random_path/dummy_non_existing_file"    

    # WHEN checking for file string
    result = get_file_status_string(str(file_not_exist))

    # THEN it should not return empty str
    assert "missing" in result[0].value_no_colors


def test_get_from_two_key():
    # GIVEN a dictionary with two keys that each have list of values
    input_dict = {"key_1":["key_1_value_1", "key_1_value_2"], "key_2": ["key_2_value_1", "key_2_value_2"]}

    # WHEN knowing the key_1_value_2 from key_1, return key_2_value_2 from key_2
    result = get_from_two_key(input_dict, from_key="key_1", by_key="key_2", by_value="key_1_value_2", default=None)
    
    # THEN retrun value should be key_2_value_2 and not None
    assert result == "key_2_value_2" 


def test_find_file_index(tmpdir):
    # GIVEN an existing bam file and its bai index file
    bam_dir = tmpdir.mkdir("temporary_path")
    bam_file = bam_dir.join("file_exists.bam")
    bam_file.write("dummy_file_content")

    bai_file = bam_dir.join("file_exists.bam.bai")
    bai_file.write("dummy_file_content")

    bai_file_2= bam_dir.join("file_exists.bai")
    bai_file_2.write("dummy_file_content")

    # WHEN finding list of bai files
    result = find_file_index(str(bam_file))
    
    # THEN return list bai file(s) as a list
    assert len(result)==2
    assert isinstance(result, list)
    assert str(bai_file) in result
    assert str(bai_file_2) in result
