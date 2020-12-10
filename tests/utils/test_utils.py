import json
import pytest
import sys
import copy
import collections

import shutil
from unittest import mock
import logging

from pathlib import Path

from BALSAMIC.utils.exc import BalsamicError, WorkflowRunError

from BALSAMIC.utils.constants import CONTAINERS_CONDA_ENV_PATH
from BALSAMIC.utils.constants import BIOINFO_TOOL_ENV
from BALSAMIC.utils.constants import REFERENCE_FILES

from BALSAMIC.utils.cli import (
    SnakeMake, CaptureStdout, iterdict, get_snakefile, createDir, write_json,
    get_config, recursive_default_dict, convert_defaultdict_to_regular_dict,
    get_file_status_string, get_from_two_key, find_file_index, merge_json,
    validate_fastq_pattern, get_panel_chrom, get_bioinfo_tools_version,
    create_fastq_symlink, get_fastq_bind_path, singularity, get_file_extension)

from BALSAMIC.utils.rule import (get_chrom, get_vcf, get_sample_type,
                                 get_picard_mrkdup, get_variant_callers,
                                 get_script_path, get_result_dir, get_threads,
                                 get_delivery_id, get_reference_output_files)

from BALSAMIC.utils.workflowscripts import get_file_contents, get_densityplot


def test_get_variant_callers_wrong_analysis_type(tumor_normal_config):
    # GIVEN a wrong analysis_type
    wrong_analysis_type = "cohort"
    workflow = "BALSAMIC"
    mutation_type = "SNV"
    mutation_class = "germline"

    # WHEN getting list of variant callers
    # THEN capture error
    with pytest.raises(WorkflowRunError):
        assert get_variant_callers(config=tumor_normal_config,
                                   analysis_type=wrong_analysis_type,
                                   workflow_solution=workflow,
                                   mutation_type=mutation_type,
                                   mutation_class=mutation_class)


def test_get_variant_callers_wrong_workflow(tumor_normal_config):
    # GIVEN a wrong workflow name
    wrong_workflow = "MIP"
    mutation_type = "SNV"
    mutation_class = "germline"
    analysis_type = "paired"

    with pytest.raises(WorkflowRunError):
        assert get_variant_callers(config=tumor_normal_config,
                                   analysis_type=analysis_type,
                                   workflow_solution=wrong_workflow,
                                   mutation_type=mutation_type,
                                   mutation_class=mutation_class)


def test_get_variant_callers_wrong_mutation_type(tumor_normal_config):
    # GIVEN a wrong workflow name
    workflow = "BALSAMIC"
    wrong_mutation_type = "INDEL"
    mutation_class = "germline"
    analysis_type = "paired"

    # WHEN getting list of variant callers
    # THEN capture error
    with pytest.raises(WorkflowRunError):
        assert get_variant_callers(config=tumor_normal_config,
                                   analysis_type=analysis_type,
                                   workflow_solution=workflow,
                                   mutation_type=wrong_mutation_type,
                                   mutation_class=mutation_class)


def test_get_variant_callers_wrong_mutation_class(tumor_normal_config):
    # GIVEN a wrong workflow name
    workflow = "BALSAMIC"
    mutation_type = "SNV"
    wrong_mutation_class = "mosaic"
    analysis_type = "paired"

    # WHEN getting list of variant callers
    # THEN capture error
    with pytest.raises(WorkflowRunError):
        assert get_variant_callers(config=tumor_normal_config,
                                   analysis_type=analysis_type,
                                   workflow_solution=workflow,
                                   mutation_type=mutation_type,
                                   mutation_class=wrong_mutation_class)


def test_get_reference_output_files():
    # GIVEN a reference genome version
    genome_ver = 'hg38'
    file_type = 'fasta'

    # WHEN getting list of valid types
    fasta_files = get_reference_output_files(REFERENCE_FILES[genome_ver],
                                             file_type)

    # THEN it should return list of file
    assert 'Homo_sapiens_assembly38.fasta' in fasta_files


def test_get_bioinfo_tools_list():
    # GIVEN a path for container path and bioinfo tool dictionary
    # WHEN getting dictionary of bioinformatic tools and their version
    bioinfo_tools_dict = get_bioinfo_tools_version(BIOINFO_TOOL_ENV, CONTAINERS_CONDA_ENV_PATH)

    # THEN assert it is a dictionary and versions are correct
    assert isinstance(bioinfo_tools_dict, dict)
    assert set(bioinfo_tools_dict["samtools"]) == set(["1.10", "1.9"])


def test_get_delivery_id():
    # GIVEN a delivery id, a dummy file string, list of tags, and a snakemake wildcard_dict
    delivery_id_to_check = '{case_name}'
    tags = ['angry_bird', 'tag_2']
    dummy_file_string = 'some_file_angry_bird_with_result.txt'
    wildcard_dict = {'case_name': 'angry_bird', "allow_missing": True}
    actual_delivery_id = 'angry_bird'

    # WHEN getting the correct delivery_id
    delivery_id_candidate = get_delivery_id(
        id_candidate=delivery_id_to_check,
        file_to_store=dummy_file_string,
        tags=tags,
        output_file_wildcards=wildcard_dict)

    # THEN correct delivery_id should be extracted
    assert delivery_id_candidate == actual_delivery_id


def test_get_file_extension_get_any_ext():
    # GIVEN a dummy file string
    dummy_file = "hassan.txt"
    actual_extension = "txt"

    # WHEN extracting the extension
    file_extension = get_file_extension(dummy_file)

    # THEN assert extension is correctly extracted
    assert file_extension == actual_extension


def test_get_file_extension_known_ext():
    # GIVEN a dummy file string with a known string
    dummy_file = "hassan.fastq.gz"
    actual_extension = "fastq.gz"

    # WHEN extracting the extension
    file_extension = get_file_extension(dummy_file)

    # THEN assert extension is correctly extracted
    assert file_extension == actual_extension


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


def test_iterdict(config_files):
    """ GIVEN a dict for iteration """
    test_dict = json.load(open(config_files['reference'], 'r'))

    # WHEN passing dict to this function
    dict_gen = iterdict(test_dict)

    # THEN it will create dict generator, we can iterate it, get the key, values as string
    for key, value in dict_gen:
        assert isinstance(key, str)
        assert isinstance(value, str)


def test_snakemake_local():
    # GIVEN required params
    snakemake_local = SnakeMake()
    snakemake_local.working_dir = "this_path/snakemake"
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
    assert "this_path/snakemake" in shell_command
    assert "--dryrun" in shell_command
    assert "--forceall" in shell_command


def test_snakemake_slurm():
    # GIVEN required params
    snakemake_slurm = SnakeMake()
    snakemake_slurm.case_name = "test_case"
    snakemake_slurm.working_dir = "this_path/snakemake"
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
    snakemake_slurm.quiet = True
    snakemake_slurm.use_singularity = True
    snakemake_slurm.singularity_bind = ["path_1", "path_2"]
    snakemake_slurm.run_analysis = True

    # WHEN calling the build command
    shell_command = snakemake_slurm.build_cmd()

    # THEN constructing snakecommand for slurm runner
    assert isinstance(shell_command, str)
    assert "worflow/variantCalling_paired" in shell_command
    assert "sample_config.json" in shell_command
    assert "this_path/snakemake" in shell_command
    assert "--dryrun" not in shell_command
    assert "sbatch.py" in shell_command
    assert "test_case" in shell_command
    assert "containers" in shell_command
    assert "--quiet" in shell_command


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
                ("single", "targeted"), ("qc", ""), ("generate_ref", ""),
                ("umi", "")]

    # WHEN asking to see snakefile for paired
    for analysis_type, sequencing_type in workflow:
        snakefile = get_snakefile(analysis_type, sequencing_type)
        pipeline = ''

        if sequencing_type in ['targeted', 'wgs', 'qc']:
            pipeline = "BALSAMIC/workflows/balsamic.smk"
        elif analysis_type == 'generate_ref':
            pipeline = "BALSAMIC/workflows/reference.smk"
        elif analysis_type == 'umi':
            pipeline = "BALSAMIC/workflows/UMIworkflow.smk"

        # THEN it should return the snakefile path
        # THEN assert file exists
        assert snakefile.startswith('/')
        assert pipeline in snakefile
        assert Path(snakefile).is_file()


def test_get_chrom(config_files):
    # Given a panel bed file
    bed_file = config_files["panel_bed_file"]
    actual_chrom = [
        '10', '11', '16', '17', '18', '19', '2', '3', '4', '6', '7', '9', 'X'
    ]

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
        ref_json = {"path": "this_path", "reference": ""}
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
    input_dict = {
        "key_1": ["key_1_value_1", "key_1_value_2"],
        "key_2": ["key_2_value_1", "key_2_value_2"]
    }

    # WHEN knowing the key_1_value_2 from key_1, return key_2_value_2 from key_2
    result = get_from_two_key(input_dict,
                              from_key="key_1",
                              by_key="key_2",
                              by_value="key_1_value_2",
                              default=None)

    # THEN retrun value should be key_2_value_2 and not None
    assert result == "key_2_value_2"


def test_find_file_index(tmpdir):
    # GIVEN an existing bam file and its bai index file
    bam_dir = tmpdir.mkdir("temporary_path")
    bam_file = bam_dir.join("file_exists.bam")
    bam_file.write("dummy_file_content")

    bai_file = bam_dir.join("file_exists.bam.bai")
    bai_file.write("dummy_file_content")

    bai_file_2 = bam_dir.join("file_exists.bai")
    bai_file_2.write("dummy_file_content")

    # WHEN finding list of bai files
    result = find_file_index(str(bam_file))

    # THEN return list bai file(s) as a list
    assert len(result) == 2
    assert isinstance(result, list)
    assert str(bai_file) in result
    assert str(bai_file_2) in result


def test_singularity_shellcmd(singularity_container):
    """test singularity shell cmd
    """

    # GIVEN a dummy command
    dummy_command = 'ls this_path'
    dummy_path_1 = 'this_path/path1'
    dummy_path_2 = 'this_path/path2'
    correct_shellcmd = 'exec --bind {} --bind {} ls this_path'.format(
        dummy_path_1, dummy_path_2)

    with mock.patch.object(shutil, 'which') as mocked:
        mocked.return_value = "/my_home/binary_path/singularity"

        # WHEN building singularity command
        shellcmd = singularity(sif_path=singularity_container,
                               cmd=dummy_command,
                               bind_paths=[dummy_path_1, dummy_path_2])

        # THEN successfully return a correct singularity cmd
        assert correct_shellcmd in shellcmd


def test_singularity_shellcmd_sif_not_exist():
    """test singularity shell cmd with non-existing file
    """

    # GIVEN a dummy command
    dummy_command = 'ls this_path'
    dummy_sif_path = '/some_path/my_sif_path_3.1415/container.sif'
    dummy_path_1 = 'this_path/path1'
    dummy_path_2 = 'this_path/path2'
    error_msg = "container file does not exist"

    # WHEN building singularity command
    # THEN successfully get error that container doesn't exist
    with mock.patch.object(shutil,
                           'which') as mocked, pytest.raises(BalsamicError,
                                                             match=error_msg):
        mocked.return_value = "/my_home/binary_path/singularity"

        singularity(sif_path=dummy_sif_path,
                    cmd=dummy_command,
                    bind_paths=[dummy_path_1, dummy_path_2])


def test_singularity_shellcmd_cmd_not_exist(singularity_container):
    """test singularity shell cmd with nonexisting singularity command
    """

    # GIVEN a dummy command
    dummy_command = 'ls this_path'
    error_msg = "singularity command does not exist"
    dummy_path_1 = 'this_path/path1'
    dummy_path_2 = 'this_path/path2'

    # WHEN building singularity command
    # THEN successfully get error if singualrity command doesn't exist
    with mock.patch.object(shutil,
                           'which') as mocked, pytest.raises(BalsamicError,
                                                             match=error_msg):
        mocked.return_value = None

        singularity(sif_path=singularity_container,
                    cmd=dummy_command,
                    bind_paths=[dummy_path_1, dummy_path_2])


def test_merge_json(config_files):
    # GIVEN a dict and json file
    ref_dict = json.load(open(config_files['reference'], 'r'))

    json_file = config_files['sample']

    # WHEN passing dict and json file to merge
    merge_dict = merge_json(ref_dict, json_file)

    # THEN It will merge both the data and return dict
    assert isinstance(merge_dict, dict)
    assert "samples" in merge_dict
    assert "reference" in merge_dict


def test_merge_json_error(config_files):
    with pytest.raises(Exception, match=r"No such file or directory"):
        # GIVEN a dict and invalid json file path
        ref_dict = json.load(open(config_files['reference'], 'r'))
        json_file = 'reference.json'

        # WHEN passing python dict and invalid json path
        # THEN it should throw OSError as FileNotFoundError
        assert merge_json(ref_dict, json_file)


def test_validate_fastq_pattern():
    #GIVEN a path to a file with correct fastq file prefix
    fastq_path_r1 = "/home/analysis/dummy_tumor_R_1.fastq.gz"
    fastq_path_r2 = "/home/analysis/dummy_normal_R_2.fastq.gz"
    #THEN it should return the correct prefix
    assert validate_fastq_pattern(fastq_path_r1) == "dummy_tumor_R"
    assert validate_fastq_pattern(fastq_path_r2) == "dummy_normal_R"

    with pytest.raises(AttributeError) as excinfo:
        #GIVEN a path to a file with incorrect fastq file prefix
        bad_fastq_path_1 = "/home/analysis/dummy_tumor.fastq.gz"
        validate_fastq_pattern(bad_fastq_path_1)
        #THEN AttributeError is raised
    assert excinfo.value

    with pytest.raises(AttributeError) as excinfo:
        #GIVEN a path to a file with incorrect fastq file prefix
        bad_fastq_path_2 = "/home/analysis/dummy_tumor_R3.fastq.gz"
        validate_fastq_pattern(bad_fastq_path_2)
        #THEN AttributeError is raised
    assert excinfo.value

    with pytest.raises(AttributeError) as excinfo:
        #GIVEN a path to a file with incorrect fastq file prefix
        bad_fastq_path_3 = "/home/analysis/dummy_tumor_R_2.bam"
        validate_fastq_pattern(bad_fastq_path_3)
        #THEN AttributeError is raised
    assert excinfo.value


def test_get_panel_chrom():
    #GIVEN a valid PANEL BED file
    panel_bed_file = 'tests/test_data/references/panel/panel.bed'
    #THEN it should return a set containing multiple unique chromosomes
    assert len(get_panel_chrom(panel_bed_file)) > 0


def test_create_fastq_symlink(tmpdir_factory, caplog):
    #GIVEN a list of valid input fastq files from test directory containing 4 files
    symlink_from_path = tmpdir_factory.mktemp("symlink_from")
    symlink_to_path = tmpdir_factory.mktemp("symlink_to")
    filenames = [
        "tumor_R_1.fastq.gz", "normal_R_1.fastq.gz", "tumor_R_2.fastq.gz",
        "normal_R_2.fastq.gz"
    ]
    successful_log = "skipping"
    casefiles = [Path(symlink_from_path, x) for x in filenames]
    for casefile in casefiles:
        casefile.touch()
    with caplog.at_level(logging.INFO):
        create_fastq_symlink(casefiles=casefiles, symlink_dir=symlink_to_path)
        #THEN destination should have 4 files
        assert len(list(Path(symlink_to_path).rglob("*.fastq.gz"))) == 4
        #THEN exception triggers log message containing "skipping"
        assert successful_log in caplog.text


def test_get_fastq_bind_path(tmpdir_factory):
    #GIVEN a list of valid input fastq filenames and test directories
    filenames = [
        "tumor_R_1.fastq.gz", "normal_R_1.fastq.gz", "tumor_R_2.fastq.gz",
        "normal_R_2.fastq.gz"
    ]
    #WHEN files are created, and symlinks are made in symlink directory
    symlink_from_path = tmpdir_factory.mktemp("symlink_from")
    symlink_to_path = tmpdir_factory.mktemp("symlink_to")
    casefiles = [Path(symlink_from_path, x) for x in filenames]
    for casefile in casefiles:
        casefile.touch()
    create_fastq_symlink(casefiles=casefiles, symlink_dir=symlink_to_path)
    #THEN function returns list containing the original parent path!
    assert get_fastq_bind_path(symlink_to_path) == [symlink_from_path]


def test_get_file_contents():
    #GIVEN a test input file
    test_file = 'tests/test_data/densityplots/dummy_file1.txt'

    # WHEN invoking function
    test_file_built = get_file_contents(test_file, 'umi')
    column_names = ['id', 'AF', 'method']

    # THEN check column names and no. of column matches
    assert all(test_file_built.columns == column_names)
    assert len(test_file_built.columns) == 3


def test_get_wrongfile_contents():
    #GIVEN a test input file
    test_wrongfile = 'tests/test_data/densityplots/dummy_wrongfile.txt'

    # WHEN invoking function
    with pytest.raises(ValueError):
        test_wrongfile_built = get_file_contents(test_wrongfile, 'umi')
        assert len(test_wrongfile_built.columns) != 3


def test_get_densityplot():
    #GIVEN prefix names, input and files
    test_file1 = 'tests/test_data/densityplots/dummy_file1.txt'
    test_file2 = 'tests/test_data/densityplots/dummy_file2.txt'
    name1 = 'testnam1'
    name2 = 'testnam2'
    out_file = 'tests/test_data/densityplots/dummy_plot.pdf'

    # WHEN invoking function out_file is created
    test_result = get_densityplot(test_file1, test_file2, name1, name2,
                                  out_file)
    test_result_name = Path(test_result).name

    # THEN check for filepaths
    assert Path(test_file1).exists()
    assert Path(test_file2).exists()
    assert Path(out_file).exists()
    assert test_result_name == "dummy_plot.pdf"
