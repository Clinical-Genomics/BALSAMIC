import json
import subprocess
import pytest
import sys
import copy
import collections

import shutil
from unittest import mock
import logging

from pathlib import Path

from _pytest.logging import LogCaptureFixture

from BALSAMIC import __version__ as balsamic_version
from BALSAMIC.utils.exc import BalsamicError, WorkflowRunError

from BALSAMIC.constants.common import CONTAINERS_CONDA_ENV_PATH, BIOINFO_TOOL_ENV
from BALSAMIC.constants.reference import REFERENCE_FILES

from BALSAMIC.utils.cli import (
    SnakeMake,
    CaptureStdout,
    iterdict,
    get_snakefile,
    createDir,
    get_config,
    get_file_status_string,
    find_file_index,
    validate_fastq_pattern,
    get_panel_chrom,
    singularity,
    get_file_extension,
    get_bioinfo_tools_version,
    convert_deliverables_tags,
    check_executable,
    job_id_dump_to_yaml,
    generate_h5,
    get_md5,
    create_md5,
    get_pon_sample_dict,
    get_input_files_path,
    get_sample_dict,
    get_fastq_info,
)
from BALSAMIC.utils.io import read_json, write_json, read_yaml

from BALSAMIC.utils.rule import (
    get_chrom,
    get_vcf,
    get_sample_type,
    get_picard_mrkdup,
    get_variant_callers,
    get_script_path,
    get_result_dir,
    get_threads,
    get_delivery_id,
    get_reference_output_files,
    get_rule_output,
    get_sample_type_from_prefix,
)



def test_get_variant_callers_wrong_analysis_type(tumor_normal_config):
    # GIVEN a wrong analysis_type
    wrong_analysis_type = "cohort"
    workflow = "BALSAMIC"
    mutation_type = "SNV"
    sequencing_type = "wgs"
    mutation_class = "germline"

    # WHEN getting list of variant callers
    # THEN capture error
    with pytest.raises(WorkflowRunError):
        assert get_variant_callers(
            config=tumor_normal_config,
            analysis_type=wrong_analysis_type,
            workflow_solution=workflow,
            mutation_type=mutation_type,
            mutation_class=mutation_class,
            sequencing_type=sequencing_type,
        )


def test_get_variant_callers_wrong_workflow(tumor_normal_config):
    # GIVEN a wrong workflow name
    wrong_workflow = "MIP"
    mutation_type = "SNV"
    mutation_class = "germline"
    sequencing_type = "wgs"
    analysis_type = "paired"

    with pytest.raises(WorkflowRunError):
        assert get_variant_callers(
            config=tumor_normal_config,
            analysis_type=analysis_type,
            workflow_solution=wrong_workflow,
            mutation_type=mutation_type,
            mutation_class=mutation_class,
            sequencing_type=sequencing_type,
        )


def test_get_variant_callers_wrong_mutation_type(tumor_normal_config):
    # GIVEN a wrong workflow name
    workflow = "BALSAMIC"
    wrong_mutation_type = "INDEL"
    mutation_class = "germline"
    sequencing_type = "wgs"
    analysis_type = "paired"

    # WHEN getting list of variant callers
    # THEN capture error
    with pytest.raises(WorkflowRunError):
        assert get_variant_callers(
            config=tumor_normal_config,
            analysis_type=analysis_type,
            workflow_solution=workflow,
            mutation_type=wrong_mutation_type,
            mutation_class=mutation_class,
            sequencing_type=sequencing_type,
        )


def test_get_variant_callers_wrong_mutation_class(tumor_normal_config):
    # GIVEN a wrong workflow name
    workflow = "BALSAMIC"
    mutation_type = "SNV"
    wrong_mutation_class = "mosaic"
    sequencing_type = "wgs"
    analysis_type = "paired"

    # WHEN getting list of variant callers
    # THEN capture error
    with pytest.raises(WorkflowRunError):
        assert get_variant_callers(
            config=tumor_normal_config,
            analysis_type=analysis_type,
            workflow_solution=workflow,
            mutation_type=mutation_type,
            mutation_class=wrong_mutation_class,
            sequencing_type=sequencing_type,
        )


def test_get_variant_callers_wrong_sequencing_type(tumor_normal_config):
    # GIVEN a wrong workflow name
    workflow = "BALSAMIC"
    mutation_type = "SNV"
    mutation_class = "somatic"
    wrong_sequencing_type = "wts"
    analysis_type = "paired"

    # WHEN getting list of variant callers
    # THEN capture error
    with pytest.raises(WorkflowRunError):
        assert get_variant_callers(
            config=tumor_normal_config,
            analysis_type=analysis_type,
            workflow_solution=workflow,
            mutation_type=mutation_type,
            mutation_class=mutation_class,
            sequencing_type=wrong_sequencing_type,
        )


def test_get_reference_output_files():
    # GIVEN a reference genome version
    genome_ver = "hg38"
    file_type = "fasta"

    # WHEN getting list of valid types
    fasta_files = get_reference_output_files(REFERENCE_FILES[genome_ver], file_type)

    # THEN it should return list of file
    assert "Homo_sapiens_assembly38.fasta" in fasta_files


def test_get_bioinfo_tools_version():
    # GIVEN a path for container path and bioinfo tool dictionary
    # WHEN getting dictionary of bioinformatic tools and their version
    bioinfo_tools_dict = get_bioinfo_tools_version(
        BIOINFO_TOOL_ENV, CONTAINERS_CONDA_ENV_PATH
    )
    observed_versions = set(bioinfo_tools_dict["samtools"])

    # THEN assert it is a dictionary and versions are correct
    assert isinstance(bioinfo_tools_dict, dict)
    assert set(observed_versions).issubset(set(["1.15.1", "1.12", "1.15", "1.9"]))


def test_get_delivery_id():
    # GIVEN a delivery id, a dummy file string, list of tags, and a snakemake wildcard_dict
    delivery_id_to_check = "{case_name}"
    tags = ["angry_bird", "tag_2"]
    dummy_file_string = "some_file_angry_bird_with_result.txt"
    wildcard_dict = {"case_name": "angry_bird", "allow_missing": True}
    actual_delivery_id = "angry_bird"

    # WHEN getting the correct delivery_id
    delivery_id_candidate = get_delivery_id(
        id_candidate=delivery_id_to_check,
        file_to_store=dummy_file_string,
        tags=tags,
        output_file_wildcards=wildcard_dict,
    )

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


def test_iterdict(reference):
    """GIVEN a dict for iteration"""
    # WHEN passing dict to this function
    dict_gen = iterdict(reference)

    # THEN it will create dict generator, we can iterate it, get the key, values as string
    for key, value in dict_gen:
        assert isinstance(key, str)
        assert isinstance(value, str)


def test_snakemake_local():
    # GIVEN required params
    snakemake_local = SnakeMake()
    snakemake_local.working_dir = "this_path/snakemake"
    snakemake_local.snakefile = "workflow/variantCalling_paired"
    snakemake_local.configfile = "sample_config.json"
    snakemake_local.run_mode = "local"
    snakemake_local.use_singularity = True
    snakemake_local.singularity_bind = ["path_1", "path_2"]
    snakemake_local.forceall = True

    # WHEN calling the build command
    shell_command = snakemake_local.build_cmd()

    # THEN it will contruct the snakemake command to run
    assert isinstance(shell_command, str)
    assert "workflow/variantCalling_paired" in shell_command
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
    snakemake_slurm.sm_opt = ("containers",)
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
        assert get_script_path(script).startswith("/")

        # THEN assert if script file exists
        assert Path(get_script_path(script)).is_file()


def test_get_snakefile():
    # GIVEN analysis_type for snakemake workflow
    workflow = [
        ("paired", "balsamic"),
        ("single", "balsamic"),
        ("generate_ref", "balsamic"),
        ("pon", "balsamic"),
        ("paired", "balsamic-qc"),
        ("single", "balsamic-qc"),
        ("paired", "balsamic-umi"),
        ("single", "balsamic-umi"),
    ]

    # WHEN asking to see snakefile for paired
    for reference_genome in ["hg19", "hg38", "canfam3"]:
        for analysis_type, analysis_workflow in workflow:
            snakefile = get_snakefile(
                analysis_type, analysis_workflow, reference_genome
            )

            pipeline = ""
            if (
                analysis_type in ["single", "paired"]
                and analysis_workflow != "balsamic-qc"
                and analysis_workflow != "balsamic-umi"
            ):
                pipeline = "BALSAMIC/workflows/balsamic.smk"
            elif analysis_type == "generate_ref" and reference_genome != "canfam3":
                pipeline = "BALSAMIC/workflows/reference.smk"
            elif analysis_type == "generate_ref" and reference_genome == "canfam3":
                pipeline = "BALSAMIC/workflows/reference-canfam3.smk"
            elif analysis_type == "pon":
                pipeline = "BALSAMIC/workflows/PON.smk"
            elif analysis_workflow == "balsamic-qc":
                pipeline = "BALSAMIC/workflows/QC.smk"

            # THEN it should return the snakefile path
            # THEN assert file exists
            assert snakefile.startswith("/")
            assert pipeline in snakefile
            assert Path(snakefile).is_file()


def test_get_chrom(config_files):
    # Given a panel bed file
    bed_file = config_files["panel_bed_file"]
    actual_chrom = [
        "10",
        "11",
        "16",
        "17",
        "18",
        "19",
        "2",
        "3",
        "4",
        "6",
        "7",
        "9",
        "X",
    ]

    # WHEN passing this bed file
    test_chrom = get_chrom(bed_file)

    # THEN It should return list of chrom presents in that bed file
    assert set(actual_chrom) == set(test_chrom)


def test_get_vcf(sample_config):
    # GIVEN a sample_config dict and a variant callers list
    variant_callers = ["tnscope", "vardict", "manta"]

    # WHEN passing args to that function
    vcf_list = get_vcf(
        sample_config, variant_callers, [sample_config["analysis"]["case_id"]]
    )

    # THEN It should return the list of vcf file names
    assert any("tnscope" in vcf_name for vcf_name in vcf_list)
    assert any("vardict" in vcf_name for vcf_name in vcf_list)
    assert any("manta" in vcf_name for vcf_name in vcf_list)


def test_get_vcf_invalid_variant_caller(sample_config):
    # GIVEN a sample_config dict and an incorrect variant callers list
    variant_callers = ["vardict", "manta", "tnhaplotyper"]

    # WHEN passing args to that function
    with pytest.raises(KeyError):
        # THEN a key error should be raised for a not supported caller
        get_vcf(sample_config, variant_callers, [sample_config["analysis"]["case_id"]])


def test_get_sample_type(sample_config):
    # GIVEN a sample_config dict, bio_type as tumor
    bio_type = "tumor"

    # WHEN calling get_sample_type with bio_type
    sample_id = get_sample_type(sample_config["samples"], bio_type)

    # THEN It should return the tumor samples id
    assert sample_id == ["S1_R"]


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
    test_stdout_message = "Message to stdout"
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
    config_file = "non_existing_config"

    # WHEN passing file names
    # THEN return the config files path
    with pytest.raises(FileNotFoundError):
        assert get_config(config_file)


def test_write_json(tmp_path, reference):
    # GIVEN a dict from sample json file
    tmp = tmp_path / "tmp"
    tmp.mkdir()
    output_json = tmp / "output.json"

    # WHEN passing dict and file name
    write_json(reference, output_json)
    output = output_json.read_text()

    # THEN It will create a json file with given dict
    for key, value in iterdict(reference):
        assert key in output
        assert value in output

    assert len(list(tmp.iterdir())) == 1


def test_write_json_error(tmp_path: Path):
    """Test JSON write error."""

    # GIVEN a dictionary to be saved in a JSON file
    ref_json = {"case": "/path/to/case", "reference": "/path/to/reference"}

    # GIVEN a directory as the output file
    with pytest.raises(Exception, match=r"Is a directory"):
        # THEN an exception should be raised
        assert write_json(ref_json, tmp_path)


def test_read_json(config_path):
    """Test data extraction from a BALSAMIC config JSON file."""

    # GIVEN a config path

    # WHEN calling the function
    config_dict = read_json(config_path)

    # THEN the config.json file should be correctly parsed
    assert type(config_dict) is dict


def test_read_json_error():
    """Test data extraction from a BALSAMIC config JSON file for an invalid path."""

    # GIVEN an incorrect config path
    config_path = "/not/a/path"

    # WHEN calling the function

    # THEN an error should raise due to an invalid file path
    try:
        read_json(config_path)
        assert False
    except FileNotFoundError as file_exc:
        assert f"The JSON file {config_path} was not found" in str(file_exc)


def test_read_yaml(metrics_yaml_path):
    """Test data extraction from a saved YAML file."""

    # GIVEN an expected output
    n_metrics = 12  # Number of expected metric

    dropout_metric = {
        "header": None,
        "id": "ACC1",
        "input": "ACC1.sorted.mrkdup.hsmetric",
        "name": "GC_DROPOUT",
        "step": "multiqc_picard_HsMetrics",
        "value": 0.027402,
        "condition": {"norm": "lt", "threshold": 1.0},
    }

    ins_size_metric = {
        "header": None,
        "id": "ACC1",
        "input": "ACC1.sorted.insertsizemetric",
        "name": "MEAN_INSERT_SIZE",
        "step": "multiqc_picard_insertSize",
        "value": 201.813054,
        "condition": None,
    }

    dups_metric = {
        "header": None,
        "id": "ACC1",
        "input": "ACC1.sorted.mrkdup.txt",
        "name": "PERCENT_DUPLICATION",
        "step": "multiqc_picard_dups",
        "value": 0.391429,
        "condition": None,
    }

    # WHEN calling the function
    requested_metrics = read_yaml(metrics_yaml_path)

    # THEN check if the data are correctly retrieved from the YAML
    assert len(requested_metrics) == n_metrics
    assert dropout_metric in requested_metrics
    assert ins_size_metric in requested_metrics
    assert dups_metric in requested_metrics


def test_read_yaml_error():
    """Test data extraction from an incorrect YAML path."""

    # GIVEN an invalid path
    yaml_path = "NOT_A_PATH"

    # THEN assert that the FileNotFoundError is raised
    try:
        read_yaml(yaml_path)
    except FileNotFoundError as file_exc:
        assert f"The YAML file {yaml_path} was not found" in str(file_exc)


def test_get_threads(config_files):
    # GIVEN cluster config file and rule name
    cluster_config = json.load(open(config_files["cluster_json"], "r"))
    rule_name = "sentieon_align_sort"

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


def test_singularity_shellcmd(balsamic_cache):
    """test singularity shell cmd"""

    # GIVEN a dummy command
    dummy_command = "ls this_path"
    dummy_path_1 = "this_path/path1"
    dummy_path_2 = "this_path/path2"
    correct_shellcmd = "exec --bind {} --bind {} ls this_path".format(
        dummy_path_1, dummy_path_2
    )
    singularity_container_sif = Path(
        balsamic_cache, balsamic_version, "containers", "align_qc", "example.sif"
    ).as_posix()

    with mock.patch.object(shutil, "which") as mocked:
        mocked.return_value = "/my_home/binary_path/singularity"

        # WHEN building singularity command
        shellcmd = singularity(
            sif_path=singularity_container_sif,
            cmd=dummy_command,
            bind_paths=[dummy_path_1, dummy_path_2],
        )

        # THEN successfully return a correct singularity cmd
        assert correct_shellcmd in shellcmd


def test_singularity_shellcmd_sif_not_exist():
    """test singularity shell cmd with non-existing file"""

    # GIVEN a dummy command
    dummy_command = "ls this_path"
    dummy_sif_path = "/some_path/my_sif_path_3.1415/container.sif"
    dummy_path_1 = "this_path/path1"
    dummy_path_2 = "this_path/path2"
    error_msg = "container file does not exist"

    # WHEN building singularity command
    # THEN successfully get error that container doesn't exist
    with mock.patch.object(shutil, "which") as mocked, pytest.raises(
        BalsamicError, match=error_msg
    ):
        mocked.return_value = "/my_home/binary_path/singularity"

        singularity(
            sif_path=dummy_sif_path,
            cmd=dummy_command,
            bind_paths=[dummy_path_1, dummy_path_2],
        )


def test_singularity_shellcmd_cmd_not_exist():
    """test singularity shell cmd with nonexisting singularity command"""

    # GIVEN a dummy command
    dummy_command = "ls this_path"
    error_msg = "singularity command does not exist"
    dummy_path_1 = "this_path/path1"
    dummy_path_2 = "this_path/path2"
    singularity_container_sif = "some_path/container.sif"

    # WHEN building singularity command
    # THEN successfully get error if singualrity command doesn't exist
    with mock.patch.object(shutil, "which") as mocked, pytest.raises(
        BalsamicError, match=error_msg
    ):
        mocked.return_value = None

        singularity(
            sif_path=singularity_container_sif,
            cmd=dummy_command,
            bind_paths=[dummy_path_1, dummy_path_2],
        )

def test_get_panel_chrom():
    # GIVEN a valid PANEL BED file
    panel_bed_file = "tests/test_data/references/panel/panel.bed"
    # THEN it should return a set containing multiple unique chromosomes
    assert len(get_panel_chrom(panel_bed_file)) > 0


def test_convert_deliverables_tags():
    """Test generation of delivery tags."""

    # GIVEN a deliverables dict and a sample config dict
    delivery_json = {
        "files": [
            {
                "path": "dummy_balsamic_run/run_tests/TN_WGS/analysis/fastq/ACC1_R_1.fp.fastq.gz",
                "path_index": [],
                "step": "fastp",
                "tag": "ACC1,read1,quality-trimmed-fastq-read1",
                "id": "ACC1",
                "format": "fastq.gz",
            },
            {
                "path": "dummy_balsamic_run/run_tests/TN_WGS/analysis/fastq/ACC1_R_2.fp.fastq.gz",
                "path_index": [],
                "step": "fastp",
                "tag": "read2,quality-trimmed-fastq-read1",
                "id": "ACC1",
                "format": "fastq.gz",
            },
            {
                "path": "dummy_balsamic_run/run_tests/TN_WGS/analysis/qc/fastp/ACC1.fastp.json",
                "path_index": [],
                "step": "fastp",
                "tag": "ACC1,json,quality-trimmed-fastq-json,tumor",
                "id": "tumor",
                "format": "json",
            },
        ]
    }
    sample_config_dict = {"samples": {"ACC1": {"type": "tumor"}}}

    # WHEN running the convert function
    delivery_json: dict = convert_deliverables_tags(
        delivery_json=delivery_json, sample_config_dict=sample_config_dict
    )

    # THEN prefix strings should be replaced with sample name
    for delivery_file in delivery_json["files"]:
        assert delivery_file["id"] == "ACC1"
        assert "ACC1" in delivery_file["tag"]
        assert "tumor" not in delivery_file["tag"]


def test_check_executable_exists():

    # GIVEN an existing executable command
    test_command = "ls"

    # WHEN calling check_executable
    # THEN it should return True
    assert check_executable(test_command)


def test_check_executable_not_existing():

    # GIVEN an existing executable command
    test_command = "twenty_twenty_was_bad"

    # WHEN calling check_executable
    # THEN it should return True
    assert not check_executable(test_command)


def test_job_id_dump_to_yaml(tmp_path):

    # GIVEN a file with one job id per line, a key (case name), and an output file name
    dummy_dir = tmp_path / "job_id_dump_dir"
    dummy_dir.mkdir()
    dummy_job_id_dump = dummy_dir / "jod_id.dump"
    dummy_job_id_dump.write_text("01234\n56789")

    dummy_name = "angrybird"

    dummy_yaml_out = dummy_dir / "jod_id.yaml"

    # WHEN creating yaml from job id dump
    job_id_dump_to_yaml(dummy_job_id_dump, dummy_yaml_out, dummy_name)

    # THEN file should exist
    assert dummy_yaml_out.exists()


def test_generate_h5(tmp_path):

    # GIVEN a job name, a path, and a job id
    dummy_path = tmp_path / "h5dir"
    dummy_path.mkdir()
    dummy_job_name = "awesome_name"
    dummy_job_id = "31415.123123"
    correct_output = Path(dummy_path, dummy_job_name + ".h5")

    # WHEN generating a h5 output
    with mock.patch.object(subprocess, "check_output") as mocked:
        actual_output = generate_h5(dummy_job_name, dummy_job_id, dummy_path)

    assert actual_output == correct_output


def test_generate_h5_capture_no_output(tmp_path):

    # GIVEN a job name, a path, and a job id
    dummy_path = tmp_path / "h5dir"
    dummy_path.mkdir()
    dummy_job_name = "awesome_name"
    dummy_job_id = "31415.123123"
    mocked_output = "sh5util: No node-step files found for jobid"
    correct_output = Path(dummy_path, dummy_job_name + ".h5")

    # WHEN generating a h5 output
    with mock.patch.object(subprocess, "check_output") as mocked:
        mocked.return_value = mocked_output.encode("utf-8")
        actual_output = generate_h5(dummy_job_name, dummy_job_id, dummy_path)

    assert actual_output == None


def test_get_md5(tmp_path):

    # GIVEN a dummy file
    dummy_dir = tmp_path / "md5"
    dummy_dir.mkdir()
    dummy_file = dummy_dir / "dummy_file.dump"
    dummy_file.write_text("Awesome Text")

    # THEN md5 returned should be
    assert get_md5(dummy_file) == "3945B39E"


def test_create_md5(tmp_path):

    # GIVEN a path to a md5 file and reference dummy files
    ref_dir = tmp_path / "references"
    ref_dir.mkdir()
    dummy_ref_file1 = ref_dir / "reference_file1.dump"
    dummy_ref_file1.write_text("Test reference1")
    dummy_ref_file2 = ref_dir / "reference_file2.dump"
    dummy_ref_file2.write_text("Test reference2")
    dummy_reference_dict = {
        "reference_dummy1": str(dummy_ref_file1),
        "reference_dummy2": str(dummy_ref_file2),
    }
    dummy_dir = tmp_path / "md5"
    dummy_dir.mkdir()
    dummy_file = dummy_dir / "dummy_file.dump"

    create_md5(dummy_reference_dict, dummy_file)

    # THEN md5 file exists
    assert dummy_file.exists()


def test_get_rule_output(snakemake_fastqc_rule):
    """Tests retrieval of existing output files from a specific workflow."""

    # GIVEN a snakemake fastqc rule object, a rule name and a list of associated wildcards
    rules = snakemake_fastqc_rule
    rule_name = "fastqc"
    output_file_wildcards = {
        "sample": ["ACC1", "tumor", "normal"],
        "case_name": "sample_tumor_only",
    }

    # THEN retrieve the output files
    output_files = get_rule_output(rules, rule_name, output_file_wildcards)

    # THEN check that the fastq files has been picked up by the function and that the tags has been correctly created
    assert len(output_files) == 2
    for file in output_files:
        # Expected file names
        assert (
            Path(file[0]).name == "concatenated_ACC1_R_1.fastq.gz"
            or Path(file[0]).name == "concatenated_ACC1_R_2.fastq.gz"
        )
        # Expected tags
        assert (
            file[3] == "1,fastqc,quality-trimmed-seq-fastqc"
            or file[3] == "2,fastqc,quality-trimmed-seq-fastqc"
        )


def test_get_sample_type_from_prefix(config_dict):
    """Test sample type extraction from a extracted config file."""

    # GIVEN a config dictionary

    # GIVEN a sample name
    sample = "ACC1"

    # WHEN calling the function
    sample_type = get_sample_type_from_prefix(config_dict, sample)

    # THEN the retrieved sample type should match the expected one
    assert sample_type == "tumor"


def test_get_sample_dict(tumor_sample_name: str, normal_sample_name: str, fastq_dir: str):
    """Tests sample dictionary retrieval."""

    try:
        samples: dict = get_sample_dict(tumor_sample_name=tumor_sample_name, normal_sample_name=normal_sample_name, fastq_path=fastq_dir)
        assert True
    except Exception:
        assert False

    assert tumor_sample_name in samples
    assert normal_sample_name in samples
    assert samples[tumor_sample_name]["type"] == "tumor"
    assert samples[normal_sample_name]["type"] == "normal"
    assert samples[tumor_sample_name]["fastq_info"]
    assert samples[normal_sample_name]["fastq_info"]

def test_get_pon_sample_dict(
    fastq_dir: str, tumor_sample_name: str, normal_sample_name: str
):
    """Tests sample PON dictionary retrieval."""

    # GIVEN a FASTQ directory

    # GIVEN the expected sample dictionary
    samples_expected: dict = {"ACC1": {"type": "normal"}, "ACC2": {"type": "normal"}}

    # WHEN retrieving PON samples
    samples: dict = get_pon_sample_dict(fastq_dir)

    # THEN the samples should be retrieved from the FASTQ directory
    assert samples == samples_expected


def test_get_input_files_path(fastq_dir: str, caplog: LogCaptureFixture):
    """Test get unlinked input files directory."""

    # GIVEN an input fast path

    # WHEN  extracting the input files common path
    input_directory: str = get_input_files_path(fastq_dir)

    # THEN the fastq directory should be returned
    assert input_directory == fastq_dir


def test_get_input_symlinked_files_path(
    fastq_dir: str, tmp_path: Path, caplog: LogCaptureFixture
):
    """Test remove symlinks from a directory."""

    # GIVEN a temporary fast path containing symlinked files
    for file in Path(fastq_dir).iterdir():
        Path(tmp_path, file.name).symlink_to(file)

    # WHEN  extracting the input files common path
    input_directory: str = get_input_files_path(str(tmp_path))

    # THEN the real fastq directory should be returned
    assert input_directory == fastq_dir
