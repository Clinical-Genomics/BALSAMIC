import pytest
import glob
import BALSAMIC
import json

from pathlib import Path
from click.testing import CliRunner
from BALSAMIC.commands.base import cli
from BALSAMIC.utils.cli import (merge_json, validate_fastq_pattern,
                                get_panel_chrom, get_bioinfo_tools_list,
                                get_sample_dict, get_sample_names,
                                create_fastq_symlink)


def test_cli(invoke_cli):
    # GIVEN I want to see version of the program
    # WHEN I am asking to see version
    result = invoke_cli(['--version'])

    # THEN It should show the version of the program
    assert BALSAMIC.__version__ in result.output


def test_config(invoke_cli):
    # GIVEN I want to see config command options
    # WHEN asking to show config options
    result = invoke_cli(['config'])

    # THEN It should show config options in result
    assert 'case' in result.output
    assert 'reference' in result.output


def test_config_case(invoke_cli):
    # GIVEN want to see config-sample params with help option
    # WHEN asking to show params for config-sample
    result = invoke_cli(['config', 'case', '--help'])

    # THEN It should show all params reuired for config-sample
    assert 'sample-id' in result.output
    assert result.exit_code == 0


def test_config_case_missing_opt(invoke_cli):
    # WHEN invoking command with missing options
    result = invoke_cli(['config', 'case'])

    # THEN It should throw missiong option error
    assert 'Error: Missing option' in result.output
    assert result.exit_code == 2


def test_report_deliver(invoke_cli):
    # WHEN invoking command with missing options
    result = invoke_cli(['report', 'deliver', '--help'])

    # THEN It should throw missiong option error
    assert result.exit_code == 0


def test_report_status(invoke_cli):
    # WHEN invoking command with missing options
    result = invoke_cli(['report', 'status', '--help'])

    # THEN It should throw missiong option error
    assert result.exit_code == 0


def test_plugins(invoke_cli):
    # GIVEN want to see config-sample params with help option
    # WHEN asking to show params for config-sample
    result = invoke_cli(['plugins', '--help'])

    # THEN It should show all params reuired for config-sample
    assert result.exit_code == 0


def test_plugins_scout(invoke_cli):
    # WHEN invoking command with missing options
    result = invoke_cli(['plugins', 'scout', '--help'])

    # THEN It should throw missiong option error
    assert result.exit_code == 0


def test_run(invoke_cli):
    # WHEN asking to options for run command
    result = invoke_cli(['run', '--help'])

    # THEN It should show all the params without any error
    assert result.exit_code == 0
    assert "analysis" in result.output


def test_run_analysis(invoke_cli):
    # WHEN invoking run analysis command
    result = invoke_cli(['run', 'analysis', '--help'])

    # THEN it should show all params without error
    assert "--analysis-type" in result.output
    assert "--snake-file" in result.output
    assert "--sample-config" in result.output
    assert "--run-mode" in result.output
    assert "--cluster-config" in result.output
    assert "--run-analysis" in result.output


def test_run_missing_opt(invoke_cli):
    # WHEN invoking run command with missing option
    result = invoke_cli(['run', 'analysis'])

    # THEN It should throw missiong option error
    assert 'Error: Missing option' in result.output
    assert result.exit_code == 2


def test_run_analysis_invalid(invoke_cli):
    # WHEN invoking run with invalid input value
    result = invoke_cli(['run', 'analysis', '--run-mode', 'foo'])

    # THEN It should throw invalid value error
    assert result.exit_code == 2
    assert 'Error: Invalid value' in result.output


def test_run_reference(invoke_cli):
    # WHEN invoking run reference command
    result = invoke_cli(['run', 'reference', '--help'])

    # THEN It should show the help message with all params
    assert "--snakefile" in result.output
    assert "--configfile" in result.output
    assert "--run-mode" in result.output
    assert "--cluster-config" in result.output
    assert "--run-analysis" in result.output
    assert result.exit_code == 0


def test_run_ref_invalid(invoke_cli):
    # WHEN invoking run reference command with invalid param
    result = invoke_cli(['run', 'reference', '--run-mode', 'foo'])

    # THEN It should throw invalid value error
    assert result.exit_code == 2
    assert 'Error: Invalid value' in result.output


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


def test_create_fastq_symlink(tmpdir_factory):
    #GIVEN a list of 2 valid input fastq files from test directory containing 4 files
    symlink_from_path = tmpdir_factory.mktemp("symlink_from")
    symlink_to_path = tmpdir_factory.mktemp("symlink_to")
    filenames = ["tumor_R_1.fastq.gz", "normal_R_1.fastq.gz", "tumor_R_2.fastq.gz", "normal_R_2.fastq.gz"]
    casefiles = [Path(symlink_from_path,x) for x in filenames]
    for casefile in casefiles:
        casefile.touch()
    #THEN destination should have 4 files
    create_fastq_symlink(casefiles=casefiles[:2], symlink_dir=symlink_to_path)
    assert len(list(Path(symlink_to_path).rglob("*.fastq.gz"))) == 4
