import json
from pathlib import Path
from BALSAMIC.utils.qc_metrics import get_qc_metrics


def test_get_metrics():

    # GIVEN analysis_path and json metric file
    analysis_path = 'tests/test_data/qc_files/analysis'
    single_json = 'multiqc_picard_HsMetrics.json'

    # WHEN invokes function
    qc_data = get_qc_metrics(analysis_path)

    # THEN assert it is as a dictionary
    assert isinstance(qc_data, dict)


def test_file_name():

    # GIVEN analysis path and json files
    analysis_path = 'tests/test_data/qc_files/analysis'
    single_json = 'multiqc_picard_HsMetrics.json'

    # WHEN looking for files
    file_name = analysis_path + "/qc/multiqc_data/" + single_json

    # THEN assert for file exits
    assert Path(file_name).exists()


def test_input_list():

    # GIVEN input list of dicts
    multiple_json = [{
        'file_name': 'test_insertsize.json',
        'required_metrics': ['coverage', 'dups']
    }, {
        'file_name': 'test_HSmetrics.json',
        'required_metrics': ['fold80', 'pct']
    }]

    # WHEN looking for dict keys
    json_keys = ['file_name', 'required_metrics']
    multiple_json_keys = []
    for j in multiple_json:
        multiple_json_keys += list(j.keys())

    # THEN assert if two lists matches
    assert list(set(multiple_json_keys)) == list(set(json_keys))


def test_json_file():

    # GIVEN
    actual_json = {
        'concatenated_ACCyyyy_XXXXXX_R': {
            'SAMPLE_NAME': 'concatenated_ACC7yyyy_XXXXXX_R',
            'FOLD_80_BASE_PENALTY': 201.0,
            'MEAN_TARGET_COVERAGE': 150.0
        }
    }

    required_metrics = ['FOLD_80_BASE_PENALTY', 'MEAN_TARGET_COVERAGE']

    # WHEN access keys in actual json
    actual_json_keys = list(actual_json[list(actual_json.keys())[0]].keys())

    # THEN assert
    assert set(required_metrics).issubset(set(actual_json_keys))


def test_sample_ids():
  
    # GIVEN sampleid
    sampleid = 'concatenated_ACCyyyy_XXXXXX_R'
    expected_sampleid = 'ACCyyyy'

    # WHEN apply split 
    built_sampleid = sampleid.split('_')[1]

    # THEN assert
    assert built_sampleid == expected_sampleid
