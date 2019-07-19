import os
import re
import json
import pytest
from datetime import datetime
from pathlib import Path

from BALSAMIC.commands.config.sample import merge_json, \
    set_panel_bed, get_output_config, get_sample_config, get_analysis_type, \
    check_exist
from BALSAMIC.commands.config.sample import configure_fastq, link_fastq
from BALSAMIC.commands.config.sample import get_fastq_path
from BALSAMIC.utils.cli import iterdict, write_json, get_config


def test_get_config(install_config):#config_files):
    # GIVEN the config files name
    files = [
        "install", "sample", "analysis_paired",
        "analysis_paired_umi", "analysis_single", "analysis_single_umi"
    ]
    # WHEN passing file names
    for file in files:
        # THEN return the config files path
        assert Path(get_config(file)).exists()


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


def test_set_panel_bed(config_files):
    # GIVEN test reference json file and panel bed file
    ref_json = json.load(open(config_files['reference'], 'r'))

    panel_bed = config_files["panel_bed_file"]

    # WHEN passing args to that function
    json_out = set_panel_bed(ref_json, panel_bed)

    # THEN It will add two more items into the dict
    assert "chrom" in json_out['panel']


def test_set_panel_bed_error(config_files):
    with pytest.raises(Exception, match=r"No such file or directory"):
        # GIVEN test reference json file and invalid panel bed file path
        ref_json = json.load(open(config_files['reference'], 'r'))
        panel_bed = "panel/panel.bed"

        # WHEN passing args to that function
        # THEN It will add two more items into the dict
        assert set_panel_bed(ref_json, panel_bed)


def test_check_exist(config_files):
    # GVIEN a file path
    reference_json = config_files['reference']

    # WHEN passing ref file path
    # THEN It will return true if the file exists in path
    assert check_exist(reference_json)


def test_check_exist_error():
    with pytest.raises(Exception, match=r"No such file or directory"):
        # GIVEN a invalid file path
        invalid_file_path = "test.txt"

        # WHEN passing a invalid file path
        # THEN It will raise the error
        assert check_exist(invalid_file_path)


def test_get_analysis_type():
    # GIVEN umi flag(boolean value) and normal fq file
    umi_true = True
    normal_valid = 'normal.fastq'
    umi_false = False
    normal_invalid = ''

    # WHEN passing values
    # THEN it will return possible analysis type
    assert get_analysis_type(normal_valid, umi_true) == 'paired_umi'
    assert get_analysis_type(normal_invalid, umi_true) == 'single_umi'
    assert get_analysis_type(normal_valid, umi_false) == 'paired'
    assert get_analysis_type(normal_invalid, umi_false) == 'single'


def test_get_output_config():
    # GIVEN a config arg and sample id
    sample_id = 'test_sample'
    config_json = 'config.json'
    _config_json = ''

    # WHEN passing values
    config = get_output_config(config_json, sample_id)
    _config = get_output_config(_config_json, sample_id)

    # THEN it will return config json file
    assert config != ''
    assert _config != ''
    assert config.split('.')[-1] == 'json'
    assert _config.split('.')[-1] == 'json'


def test_get_sample_config(config_files):
    # GIVEN a sample config json with sample id, analysis dir, analysis type
    sample_id = 'sample1'
    sample_config = config_files['sample']
    analysis_dir = './'
    analysis_type = 'paired'

    # WHEN passing args into the function
    sample_config = get_sample_config(sample_config, sample_id, analysis_dir,
                                      analysis_type)
    date = datetime.now().strftime("%Y-%m-%d %H:%M")

    # THEN it will return the updated python dict with given values
    assert isinstance(sample_config, dict)
    assert sample_id in sample_config['analysis'].values()
    assert date in sample_config['analysis'].values()


def test_configure_fastq(sample_config, tmp_path):
    # GIVEN sample_config, normal and tumor fastq files
    fastq_src = os.path.join(sample_config['analysis']['analysis_dir'], 'fastq')
    tumor = os.path.join(fastq_src, 'S1_R_1.fastq.gz')
    fastq_prefix = ''
    fastq_dir = tmp_path / "output"
    fastq_dir.mkdir()

    # WHEN invoking configure fastq with required params
    tumor_str = configure_fastq(fastq_dir, tumor, fastq_prefix)

    # THEN It should return sample prefixes
    assert tumor_str in sample_config['samples']
    assert len(list(tmp_path.iterdir())) == 1


def test_link_fastq(sample_config, tmp_path):
    # GIVEN list of fastq files and destination to create symlink
    fastq_src = os.path.join(sample_config['analysis']['analysis_dir'], 'fastq')
    fastq_files = [os.path.join(fastq_src, file) for file in os.listdir(fastq_src)]
    dest_dir = tmp_path / "output"
    dest_dir.mkdir()

    # WHEN calling link fastq function
    link_fastq(fastq_files, dest_dir)

    # THEN it should create symlinks for fastq files
    assert len(list(dest_dir.iterdir())) == 4


def test_link_fastq_error(sample_config, tmp_path):
    # GIVEN a invalid fastq files list
    with pytest.raises(Exception) as e:
        fastq_src = os.path.join(sample_config['analysis']['analysis_dir'], 'fastq')
        fastq_files = [os.path.join(fastq_src, file) for file in os.listdir(fastq_src)]
        dest_dir = tmp_path / "output"
        dest_dir.mkdir()

        # WHEN calling link fastq
        link_fastq(fastq_files, dest_dir)

        # THEN It should return exception
        link_fastq(fastq_files, dest_dir)

        assert e.value


def test_get_fastq_path(sample_fastq):
    # GIVEN a fastq file path, and file pattern
    fastq_prefix = ''
    fastq_file = sample_fastq['fastq_valid']
    fq_pattern = re.compile(r"R_[12]" + fastq_prefix + ".fastq.gz$")

    # WHEN invoking get fastq path
    fq_str, fq_path = get_fastq_path(fastq_file, fq_pattern)

    # THEN It should return fileprefix and fastq_path
    assert fq_str == 'S1_R'
    assert os.path.dirname(fastq_file) in fq_path


def test_get_fastq_path_error(sample_fastq):
    with pytest.raises(Exception) as error:
        # GIVEN a fastq file path, and file pattern
        fastq_prefix = ''
        fastq_file = 'S1_R_1.fastq.gz'
        fq_pattern = re.compile(r"R_[12]" + fastq_prefix + ".fastq.gz$")

        # WHEN invoking get fastq path
        # THEN It should return fileprefix and fastq_path
        assert get_fastq_path(fastq_file, fq_pattern)
        assert 'FileNotFoundError' in error.value


def test_get_fqpath_mismatch_error(sample_fastq):
    with pytest.raises(Exception) as error:
        # GIVEN a fastq file path, and file pattern
        fastq_prefix = ''
        fastq_file = sample_fastq['fastq_invalid']
        fq_pattern = re.compile(r"R_[12]" + fastq_prefix + ".fastq.gz$")

        # WHEN invoking get fastq path
        # THEN It should return fileprefix and fastq_path
        assert get_fastq_path(fastq_file, fq_pattern)
        assert 'AttributeError' in error.value


