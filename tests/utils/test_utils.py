"""Test helper functions."""

import json
import logging
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List
from unittest import mock

import click
import pytest
import yaml
from _pytest.logging import LogCaptureFixture
from _pytest.tmpdir import TempPathFactory

from BALSAMIC.commands.config.case import case_config
from BALSAMIC.constants.analysis import BIOINFO_TOOL_ENV, SampleType, SequencingType
from BALSAMIC.constants.cache import CacheVersion
from BALSAMIC.constants.cluster import ClusterConfigType
from BALSAMIC.constants.constants import FileType
from BALSAMIC.constants.paths import CONTAINERS_DIR
from BALSAMIC.models.config import ConfigModel, FastqInfoModel, SampleInstanceModel
from BALSAMIC.utils.cli import (
    CaptureStdout,
    check_executable,
    convert_deliverables_tags,
    createDir,
    find_file_index,
    generate_h5,
    get_analysis_fastq_files_directory,
    get_bioinfo_tools_version,
    get_config_path,
    get_fastq_info,
    get_file_extension,
    get_file_status_string,
    get_panel_chrom,
    get_pon_sample_list,
    get_resolved_fastq_files_directory,
    get_sample_list,
    get_snakefile,
    validate_cache_version,
    validate_exome_option,
    validate_umi_min_reads,
)
from BALSAMIC.utils.exc import BalsamicError, WorkflowRunError
from BALSAMIC.utils.io import (
    read_json,
    read_vcf_file,
    read_yaml,
    write_finish_file,
    write_json,
    write_sacct_to_yaml,
    write_yaml,
)
from BALSAMIC.utils.rule import (
    get_delivery_id,
    get_fastp_parameters,
    get_result_dir,
    get_rule_output,
    get_sample_type_from_sample_name,
    get_script_path,
    get_threads,
    get_variant_callers,
    get_vcf,
)
from BALSAMIC.utils.utils import (
    get_absolute_paths_dict,
    get_relative_paths_dict,
    remove_unnecessary_spaces,
)


def test_remove_unnecessary_spaces():
    """Tests removal of unnecessary spaces from a string."""

    # GIVEN a string with unnecessary spaces
    string: str = "  Developing Balsamic   brings  me    joy "

    # WHEN calling the function
    formatted_string: str = remove_unnecessary_spaces(string)

    # THEN the extra spaces are removed
    assert formatted_string == "Developing Balsamic brings me joy"


def test_relative_paths_dict(session_tmp_path: Path, reference_file: Path):
    """Test return of a dictionary with relative paths to a provided base path."""

    # GIVEN a base path and a dictionary with absolute paths
    absolute_paths: Dict[str, Path] = {"reference": reference_file}

    # WHEN generating the relative paths dictionary
    relative_paths: Dict[str, str] = get_relative_paths_dict(
        base_path=session_tmp_path, data=absolute_paths
    )

    # THEN a dictionary with relative paths should be returned
    assert relative_paths == {"reference": reference_file.name}


def test_absolute_paths_dict(session_tmp_path: Path, reference_file: Path):
    """Test return of a dictionary with relative paths to a provided base path."""

    # GIVEN a base path and a dictionary with relative paths
    relative_paths: Dict[str, Path] = {"reference": Path(reference_file.name)}

    # WHEN generating the resolved paths dictionary
    absolute_paths: Dict[str, Path] = get_absolute_paths_dict(
        base_path=session_tmp_path, data=relative_paths
    )

    # THEN a dictionary with relative paths should be returned
    assert absolute_paths == {"reference": reference_file}


def test_get_pon_sample_dict(pon_config_dict_w_fastq: Dict):
    """Tests sample PON dictionary retrieval."""

    # GIVEN a FASTQ directory
    fastq_dir_pon = pon_config_dict_w_fastq["analysis"]["fastq_path"]
    # WHEN retrieving PON samples
    samples: List[SampleInstanceModel] = get_pon_sample_list(fastq_dir_pon)

    # THEN the samples should be retrieved from the FASTQ directory
    # And match the expected structure of pre-designed PON-config
    for sample_dict in samples:
        assert sample_dict in pon_config_dict_w_fastq["samples"]


def test_get_variant_callers_wrong_analysis_type(tumor_normal_config: Dict):
    # GIVEN a wrong analysis_type
    wrong_analysis_type = "cohort"
    workflow = "BALSAMIC"
    mutation_type = "SNV"
    sequencing_type = "wgs"
    mutation_class = "germline"

    # WHEN getting list of variant callers
    with pytest.raises(WorkflowRunError):
        # THEN capture error
        get_variant_callers(
            config=tumor_normal_config,
            analysis_type=wrong_analysis_type,
            workflow_solution=workflow,
            mutation_type=mutation_type,
            mutation_class=mutation_class,
            sequencing_type=sequencing_type,
        )


def test_get_variant_callers_wrong_workflow(tumor_normal_config: Dict):
    # GIVEN a wrong workflow name
    wrong_workflow = "MIP"
    mutation_type = "SNV"
    mutation_class = "germline"
    sequencing_type = "wgs"
    analysis_type = "paired"

    # WHEN getting list of variant callers
    with pytest.raises(WorkflowRunError):
        # THEN capture error
        get_variant_callers(
            config=tumor_normal_config,
            analysis_type=analysis_type,
            workflow_solution=wrong_workflow,
            mutation_type=mutation_type,
            mutation_class=mutation_class,
            sequencing_type=sequencing_type,
        )


def test_get_variant_callers_wrong_mutation_type(tumor_normal_config: Dict):
    # GIVEN a wrong workflow name
    workflow = "BALSAMIC"
    wrong_mutation_type = "INDEL"
    mutation_class = "germline"
    sequencing_type = "wgs"
    analysis_type = "paired"

    # WHEN getting list of variant callers
    with pytest.raises(WorkflowRunError):
        # THEN capture error
        get_variant_callers(
            config=tumor_normal_config,
            analysis_type=analysis_type,
            workflow_solution=workflow,
            mutation_type=wrong_mutation_type,
            mutation_class=mutation_class,
            sequencing_type=sequencing_type,
        )


def test_get_variant_callers_wrong_mutation_class(tumor_normal_config: Dict):
    # GIVEN a wrong workflow name
    workflow = "BALSAMIC"
    mutation_type = "SNV"
    wrong_mutation_class = "mosaic"
    sequencing_type = "wgs"
    analysis_type = "paired"

    # WHEN getting list of variant callers
    with pytest.raises(WorkflowRunError):
        # THEN capture error
        get_variant_callers(
            config=tumor_normal_config,
            analysis_type=analysis_type,
            workflow_solution=workflow,
            mutation_type=mutation_type,
            mutation_class=wrong_mutation_class,
            sequencing_type=sequencing_type,
        )


def test_get_variant_callers_wrong_sequencing_type(tumor_normal_config: Dict):
    # GIVEN a wrong workflow name
    workflow = "BALSAMIC"
    mutation_type = "SNV"
    mutation_class = "somatic"
    wrong_sequencing_type = "wts"
    analysis_type = "paired"

    # WHEN getting list of variant callers
    with pytest.raises(WorkflowRunError):
        # THEN capture error
        get_variant_callers(
            config=tumor_normal_config,
            analysis_type=analysis_type,
            workflow_solution=workflow,
            mutation_type=mutation_type,
            mutation_class=mutation_class,
            sequencing_type=wrong_sequencing_type,
        )


def test_get_bioinfo_tools_version():
    """Test bioinformatics tools and version extraction."""

    # GIVEN a tools dictionary
    bioinfo_tools: dict = get_bioinfo_tools_version(BIOINFO_TOOL_ENV, CONTAINERS_DIR)

    # THEN assert that the versions are correctly retrieved
    assert set(bioinfo_tools["picard"]).issubset({"2.27.1"})
    assert set(bioinfo_tools["samtools"]).issubset({"1.15.1", "1.12", "1.15", "1.9"})


def test_get_bioinfo_pip_tools_version():
    """Test bioinformatics tools and version extraction for a PIP specific tool."""

    # GIVEN a tools dictionary
    bioinfo_tools: dict = get_bioinfo_tools_version(BIOINFO_TOOL_ENV, CONTAINERS_DIR)

    # THEN assert that the PIP specific packages are correctly retrieved
    assert set(bioinfo_tools["cnvkit"]).issubset({"0.9.10"})


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
    dummy_file = f"dummy.{FileType.FASTQ}.{FileType.GZ}"
    actual_extension = f"{FileType.FASTQ}.{FileType.GZ}"

    # WHEN extracting the extension
    file_extension = get_file_extension(dummy_file)

    # THEN assert extension is correctly extracted
    assert file_extension == actual_extension


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
    for _reference_genome in ["hg19", "hg38", "canfam3"]:
        for analysis_type, analysis_workflow in workflow:
            snakefile = get_snakefile(analysis_type, analysis_workflow)

            pipeline = ""
            if (
                analysis_type in ["single", "paired"]
                and analysis_workflow != "balsamic-qc"
                and analysis_workflow != "balsamic-umi"
            ):
                pipeline = "BALSAMIC/workflows/balsamic.smk"
            elif analysis_type == "generate_ref":
                pipeline = "BALSAMIC/workflows/reference.smk"
            elif analysis_type == "pon":
                pipeline = "BALSAMIC/workflows/PON.smk"
            elif analysis_workflow == "balsamic-qc":
                pipeline = "BALSAMIC/workflows/QC.smk"

            # THEN it should return the snakefile path
            # THEN assert file exists
            assert snakefile.startswith("/")
            assert pipeline in snakefile
            assert Path(snakefile).is_file()


def test_get_vcf(sample_config: Dict):
    # GIVEN a sample_config dict and a variant callers list
    variant_callers = ["tnscope", "manta"]

    # WHEN passing args to that function
    vcf_list = get_vcf(
        sample_config, variant_callers, [sample_config["analysis"]["case_id"]]
    )

    # THEN It should return the list of vcf file names
    assert any("tnscope" in vcf_name for vcf_name in vcf_list)
    assert any("manta" in vcf_name for vcf_name in vcf_list)


def test_get_vcf_invalid_variant_caller(sample_config):
    # GIVEN a sample_config dict and an incorrect variant callers list
    variant_callers = ["manta", "tnhaplotyper"]

    # WHEN passing args to that function
    with pytest.raises(KeyError):
        # THEN a key error should be raised for a not supported caller
        get_vcf(sample_config, variant_callers, [sample_config["analysis"]["case_id"]])


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


def test_get_result_dir(sample_config: Dict):
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


def test_get_config_path(cluster_analysis_config_path: str):
    """Test return of a config path given its type."""

    # GIVEN an analysis config path

    # WHEN retrieving the cluster analysis configuration
    cluster_analysis: Path = get_config_path(ClusterConfigType.ANALYSIS)

    # THEN an analysis cluster json should be returned
    assert cluster_analysis.exists()
    assert cluster_analysis.as_posix() == cluster_analysis_config_path


def test_write_json(tmp_path, reference):
    # GIVEN a dict from sample json file
    tmp = tmp_path / "tmp"
    tmp.mkdir()
    output_json = tmp / "output.json"

    # WHEN passing dict and file name
    write_json(reference, output_json)
    output = output_json.read_text()

    # THEN It will create a json file with given dict
    for key, value in reference.items():
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
        write_json(ref_json, tmp_path)


def test_read_json(config_path: str):
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
        "input": "ACC1.dedup.realign.hsmetric",
        "name": "GC_DROPOUT",
        "step": "multiqc_picard_HsMetrics",
        "value": 0.027402,
        "condition": {"norm": "lt", "threshold": 1.0},
    }

    ins_size_metric = {
        "header": None,
        "id": "ACC1",
        "input": "ACC1.dedup.realign.insertsizemetric",
        "name": "MEAN_INSERT_SIZE",
        "step": "multiqc_picard_insertSize",
        "value": 201.813054,
        "condition": None,
    }

    dups_metric = {
        "header": None,
        "id": "ACC1",
        "input": "tumor.ACC1.dedup.metrics",
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


def test_write_yaml(metrics_yaml_path: str, tmp_path: Path):
    """Tests write yaml file."""

    # GIVEN a yaml file

    # GIVEN a file path to write to
    yaml_file: Path = Path(tmp_path, "write_yaml.yaml")

    # WHEN reading the yaml file
    metrics_data: Dict[str, Any] = read_yaml(metrics_yaml_path)

    # WHEN writing the yaml file from dict
    write_yaml(data=metrics_data, file_path=yaml_file.as_posix())

    # THEN assert that a file was successfully created
    assert Path.exists(yaml_file)

    # WHEN reading it as a yaml
    written_metrics_data: dict = read_yaml(yaml_file.as_posix())

    # THEN assert that all data is kept
    assert written_metrics_data == metrics_data


def test_get_threads(cluster_analysis_config_path: str):
    # GIVEN cluster config file and rule name
    cluster_config = json.load(open(cluster_analysis_config_path, "r"))
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


def test_get_panel_chrom():
    # GIVEN a valid PANEL BED file
    panel_bed_file = "tests/test_data/references/panel/panel.bed"
    # THEN it should return a set containing multiple unique chromosomes
    assert len(get_panel_chrom(panel_bed_file)) > 0


def test_convert_deliverables_tags(tumor_normal_fastq_info_correct: List[Dict]):
    """Test generation of delivery tags."""

    # GIVEN a deliverables dict and a sample config dict
    delivery_json = [
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

    sample_config_dict = {"samples": tumor_normal_fastq_info_correct}

    # WHEN running the convert function
    delivery_json: dict = convert_deliverables_tags(
        delivery_json=delivery_json, sample_config_dict=sample_config_dict
    )

    # THEN prefix strings should be replaced with sample name
    for delivery_file in delivery_json:
        assert "ACC1" in delivery_file["tag"]
        assert "tumor" not in delivery_file["tag"]
        assert delivery_file["id"] == "ACC1"


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


def test_get_sample_type_from_sample_name(config_dict: Dict):
    """Test sample type extraction from a extracted config file."""

    # GIVEN a config dictionary

    # GIVEN a sample name
    sample_name = "ACC1"

    # WHEN calling the function
    sample_type = get_sample_type_from_sample_name(config_dict, sample_name)

    # THEN the retrieved sample type should match the expected one
    assert sample_type == SampleType.TUMOR


def test_get_rule_output(snakemake_bcftools_filter_tnscope_research_tumor_only):
    """Tests retrieval of existing output files from a specific workflow."""

    # GIVEN a snakemake fastqc rule object, a rule name and a list of associated wildcards
    rules = snakemake_bcftools_filter_tnscope_research_tumor_only
    rule_name = "bcftools_filter_tnscope_research_tumor_only"
    output_file_wildcards = {
        "sample": ["ACC1", "tumor", "normal"],
        "case_name": "sample_tumor_only",
        "var_type": ["CNV", "SNV", "SV"],
    }
    # THEN retrieve the output files
    output_files = get_rule_output(rules, rule_name, output_file_wildcards)

    # THEN check that the output files and tags are correctly retrieved
    assert len(output_files) == 2
    for file in output_files:
        # Expected file names
        assert (
            Path(file[0]).name
            == "SNV.somatic.sample_tumor_only.tnscope.research.filtered.pass.vcf.gz"
        )
        # Expected tags
        assert (
            file[3]
            == "SNV,sample-tumor-only,vcf-pass-tnscope,research-vcf-pass-tnscope"
        )


def test_get_resolved_fastq_files_directory(fastq_dir: str):
    """Test get fastq directory for unlinked fastqs."""

    # GIVEN an input fastq path

    # WHEN  extracting the input files common path
    input_dir: str = get_resolved_fastq_files_directory(fastq_dir)

    # THEN the fastq directory should be returned
    assert input_dir == fastq_dir


def test_get_resolved_fastq_files_directory_symlinked_files(
    fastq_dir: str, tmp_path: Path
):
    """Test get fastq directory for symlinked files."""

    # GIVEN a temporary fastq path containing symlinked files
    for file in Path(fastq_dir).iterdir():
        Path(tmp_path, file.name).symlink_to(file)

    # WHEN  extracting the input files common path
    input_dir: str = get_resolved_fastq_files_directory(str(tmp_path))

    # THEN the real fastq directory should be returned
    assert input_dir == fastq_dir


def test_write_finish_file(json_file: Path):
    """Test finish analysis completion file generation."""

    # GIVEN a file path to write to

    # WHEN writing a json file after an analysis has been completed
    write_finish_file(file_path=json_file.as_posix())

    # THEN assert that a file was successfully created
    assert Path.exists(json_file)


def test_get_analysis_fastq_files_directory(fastq_dir: str):
    """Test get analysis fastq directory when it already exists in case folder."""

    # GIVEN an input fastq path

    # WHEN getting the analysis fastq directory
    input_dir: str = get_analysis_fastq_files_directory(
        case_dir=Path(fastq_dir).parents[1].as_posix(), fastq_path=fastq_dir
    )

    # THEN the original fastq directory should be returned
    assert input_dir == fastq_dir


def test_get_analysis_fastq_files_directory_exception(
    fastq_dir: str,
    case_id_tumor_only: str,
    tmp_path_factory: TempPathFactory,
    caplog: LogCaptureFixture,
):
    """Test get analysis fastq directory when it already exists in case folder but another path is provided."""
    caplog.set_level(logging.INFO)

    # GIVEN an input fastq path and an external case directory
    case_dir: str = tmp_path_factory.mktemp(case_id_tumor_only).as_posix()

    # WHEN getting the analysis fastq directory twice
    _input_dir: str = get_analysis_fastq_files_directory(
        case_dir=case_dir, fastq_path=fastq_dir
    )
    input_dir: str = get_analysis_fastq_files_directory(
        case_dir=case_dir, fastq_path=fastq_dir
    )

    # THEN the fastq directory should be located inside the case directory and the linking should have been skipped
    assert input_dir == Path(case_dir, "fastq").as_posix()
    assert "Skipping linking" in caplog.text


def test_get_analysis_fastq_files_directory_no_fastqs(
    fastq_dir: str, tmp_path_factory: TempPathFactory, case_id_tumor_only: str
):
    """Test get analysis fastq directory when the provided fastq directory is outside the case folder."""

    # GIVEN an external input fastq path and a case directory
    case_dir: str = tmp_path_factory.mktemp(case_id_tumor_only).as_posix()

    # WHEN getting the analysis fastq directory
    input_dir: str = get_analysis_fastq_files_directory(
        case_dir=case_dir, fastq_path=fastq_dir
    )

    # THEN the fastq directory should be located inside the case directory
    assert input_dir == Path(case_dir, "fastq").as_posix()

    # THEN the case fast files should have been linked to the provided fastq directory
    for fastq in Path(input_dir).iterdir():
        assert fastq.is_symlink()
        assert fastq.resolve().is_file()
        assert fastq_dir == fastq.resolve().parent.as_posix()


def test_get_sample_list(
    tumor_sample_name: str, normal_sample_name: str, fastq_dir_tumor_normal: str
):
    """Tests sample dictionary retrieval."""

    samples: List = get_sample_list(
        tumor_sample_name=tumor_sample_name,
        normal_sample_name=normal_sample_name,
        fastq_path=fastq_dir_tumor_normal,
    )
    assert samples[0]["name"] == tumor_sample_name
    assert samples[1]["name"] == normal_sample_name
    assert samples[0]["type"] == "tumor"
    assert samples[1]["type"] == "normal"
    assert samples[0]["fastq_info"]
    assert samples[1]["fastq_info"]


def test_get_fastq_info(tumor_sample_name: str, fastq_dir_tumor_only: str):
    """Validates that get_fastq_info assigns fastq info as expected."""
    # GIVEN a fastq_dir and sample name ACC1

    # WHEN calling the get_fastq_info function
    fastq_dict = get_fastq_info(tumor_sample_name, fastq_dir_tumor_only)

    fwd1_expected = f"{fastq_dir_tumor_only}/1_171015_HJ7TLDSX5_ACC1_XXXXXX_1.fastq.gz"
    rev1_expected = f"{fastq_dir_tumor_only}/1_171015_HJ7TLDSX5_ACC1_XXXXXX_2.fastq.gz"
    fwd2_expected = f"{fastq_dir_tumor_only}/2_171015_HJ7TLDSX5_ACC1_XXXXXX_1.fastq.gz"
    rev2_expected = f"{fastq_dir_tumor_only}/2_171015_HJ7TLDSX5_ACC1_XXXXXX_2.fastq.gz"
    fastq_info1_expected = FastqInfoModel(fwd=fwd1_expected, rev=rev1_expected)
    fastq_info2_expected = FastqInfoModel(fwd=fwd2_expected, rev=rev2_expected)

    # THEN check that the fastq_dict matches the expected fastq_dict
    expected_fastq_dict = {
        "1_171015_HJ7TLDSX5_ACC1_XXXXXX": fastq_info1_expected.model_dump(
            exclude_none=True
        ),
        "2_171015_HJ7TLDSX5_ACC1_XXXXXX": fastq_info2_expected.model_dump(
            exclude_none=True
        ),
    }
    assert fastq_dict == expected_fastq_dict


def test_get_fastq_info_symlink(tumor_sample_name: str, fastq_dir_symlinked: str):
    """Test symlinked fast info included in samples dictionary."""

    # GIVEN a fastq_dir and sample name

    # WHEN calling the get_fastq_info function
    fastq_dict: Dict[str, FastqInfoModel] = get_fastq_info(
        tumor_sample_name, fastq_dir_symlinked
    )

    # THEN verify that the resolved file links have been included in the samples dictionary
    assert "fwd_resolved" in next(iter(fastq_dict.values()))
    assert "rev_resolved" in next(iter(fastq_dict.values()))


def test_get_fastq_info_empty_fastq_dir(tumor_sample_name: str, empty_fastq_dir: str):
    """Tests if get_fastq_info correctly reports errors of not finding any fastq-files"""
    # GIVEN an empty fastq_dir and a sample name

    # WHEN calling get_fastq_info
    # THEN the following error should be found
    with pytest.raises(
        BalsamicError, match=f"No fastqs found for: {tumor_sample_name}"
    ):
        get_fastq_info(tumor_sample_name, empty_fastq_dir)


def test_get_fastq_info_double_assigned_fastq_pattern(
    tumor_sample_name: str, fastq_dir_tumor_duplicate_fastq_patterns: str
):
    """Tests if get_fastq_info correctly reports error of finding double-assigned fastq-patterns"""
    # GIVEN an empty fastq_dir and a sample name

    # WHEN calling get_fastq_info
    # THEN the following error should be found
    with pytest.raises(BalsamicError, match="Fastq name conflict. Fastq pair pattern"):
        get_fastq_info(tumor_sample_name, fastq_dir_tumor_duplicate_fastq_patterns)


def test_get_fastp_parameters(balsamic_model: ConfigModel):
    """Validate correct retrieval of WGS and TGA specific fastp parameters."""

    # GIVEN WGS config
    balsamic_model.analysis.sequencing_type = SequencingType.WGS
    fastp_params_wgs = get_fastp_parameters(balsamic_model)
    # THEN adapter trimming should be active in qual trim params
    assert "--disable_adapter_trimming" not in fastp_params_wgs["fastp_trim_qual"]
    # THEN quality trimming should be active in adapter trim params
    assert "--disable_quality_filtering" not in fastp_params_wgs["fastp_trim_adapter"]

    # GIVEN TGA config
    balsamic_model.analysis.sequencing_type = SequencingType.TARGETED
    fastp_params_tga = get_fastp_parameters(balsamic_model)
    # THEN adapter trimming should NOT be active in qual trim params
    assert "--disable_adapter_trimming" in fastp_params_tga["fastp_trim_qual"]
    # THEN quality trimming should NOT be active in adapter trim params
    assert "--disable_quality_filtering" in fastp_params_tga["fastp_trim_adapter"]


def test_validate_exome_option(panel_bed_file: str):
    # GIVEN that a panel bedfile has been supplied and exome parameter set to true
    ctx = click.Context(case_config)
    ctx.params["panel_bed"] = panel_bed_file
    # WHEN validating exome option
    # THEN exome argument should be correctly set
    assert validate_exome_option(ctx, click.Parameter, True) == True

    # GIVEN that a panel bedfile has NOT been supplied and exome parameter set to true
    ctx.params["panel_bed"] = None

    # WHEN validating exome option
    # THEN a bad parameter error should be raised
    with pytest.raises(click.BadParameter):
        validate_exome_option(ctx, click.Parameter, True)


def test_validate_cache_version_develop():
    """Test develop cache version validation."""

    # GIVEN a develop cache version
    cli_version: str = CacheVersion.DEVELOP

    # WHEN validating the provided version
    version: str = validate_cache_version(click.Context, click.Parameter, cli_version)

    # THEN the correct version should be returned
    assert version == CacheVersion.DEVELOP


def test_validate_cache_version_release():
    """Test release cache version validation."""

    # GIVEN a release cache version
    cli_version: str = "1.2.3"

    # WHEN validating the provided version
    version: str = validate_cache_version(click.Context, click.Parameter, cli_version)

    # THEN the correct version should be returned
    assert version == cli_version


def test_validate_cache_version_non_digit():
    """Test non digit release cache version validation."""

    # GIVEN an incorrect release cache version
    cli_version: str = "a.b.c"

    # WHEN validating the provided version

    # THEN a bad parameter error should be raised
    with pytest.raises(click.BadParameter):
        validate_cache_version(click.Context, click.Parameter, cli_version)


def test_validate_cache_version_wrong_format():
    """Test wrong format release cache version validation."""

    # GIVEN an incorrect release cache version
    cli_version: str = "1.2"

    # WHEN validating the provided version

    # THEN a bad parameter error should be raised
    with pytest.raises(click.BadParameter):
        validate_cache_version(click.Context, click.Parameter, cli_version)


@pytest.mark.parametrize(
    "umi_min_reads, should_fail",
    [("3,2,1", False), ("1,0,0", False), ("3,2", True), ("a,b,c", True), (None, False)],
)
def test_validate_umi_min_reads(umi_min_reads: str, should_fail: bool):
    """Test UMI min reads option validation."""
    if should_fail:
        with pytest.raises(click.BadParameter):
            validate_umi_min_reads(
                _ctx=click.Context, _param=click.Parameter, umi_min_reads=umi_min_reads
            )
    else:
        validate_umi_min_reads(
            _ctx=click.Context, _param=click.Parameter, umi_min_reads=umi_min_reads
        )


def test_read_vcf(vcf_file_path, vcf_file_gz_path):
    """Test correct reading of VCF files."""
    # GIVEN a path to two identical VCF files, one gzipped and one not

    # WHEN reading VCF file into list of strings
    vcf_contents: List[str] = read_vcf_file(vcf_file_path)
    vcf_contents_from_gz: List[str] = read_vcf_file(vcf_file_gz_path)

    # THEN first and last line of the list should match the expected values
    first_line = vcf_contents[0]
    last_line = vcf_contents[-1]
    first_line_expect = "##fileformat=VCFv4.2"
    last_line_expect = "20\t13417\t.\tC\tCGAGA\t-0.00\tLowQual\tAC=0;AF=0;AN=2;DP=46;ExcessHet=3.0103;FS=0.000;MLEAC=0;MLEAF=0;MQ=30.92;SOR=0.307\tGT:AD:DP:GQ:PL\t0/0:46,0:46:99:0,139,2084"

    assert first_line == first_line_expect
    assert last_line == last_line_expect

    # THEN the contents from the gzipped VCF should match the non-gzipped VCF
    assert vcf_contents == vcf_contents_from_gz


def test_write_sacct_to_yaml(case_id_tumor_only: str, sacct_file: Path, tmp_path: Path):
    """Tests parsing of the sacct file to create a job IDs yaml file."""

    # GIVEN a sacct file and an output yaml file
    yaml_file_path: Path = Path(tmp_path, "slurm_jobids.yaml")

    # WHEN writing the yaml file
    write_sacct_to_yaml(
        case_id=case_id_tumor_only,
        sacct_file_path=sacct_file,
        yaml_file_path=yaml_file_path,
    )

    # THEN the content should be a list of job IDs associated with a case_id
    assert yaml_file_path.is_file()
    with open(yaml_file_path, "r") as file:
        data: Dict[str, List[str]] = yaml.safe_load(file)
    assert case_id_tumor_only in data
