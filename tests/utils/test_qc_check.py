import json

from pathlib import Path

from BALSAMIC.utils.cli import write_json
from BALSAMIC.utils.qc_check import read_json, read_qc_table, check_qc_criteria, write_output
from BALSAMIC.utils.qc_check import get_bait_name, get_sample_name, get_qc_criteria, failed_qc
from BALSAMIC.utils.qc_check import get_qc_value, check_qc_criteria_simple
from BALSAMIC.utils.constants import HSMETRICS_QC_CHECK


def test_read_json():
    # GIVEN the file exist
    hs_metrics_path = "tests/test_data/qc_files/multiqc_picard_HsMetrics.json"

    # WHEN reading the file
    df = read_json(hs_metrics_path)

    # THEN check if the file contains any values
    bol_list = df.any().tolist()
    for n in range(len(bol_list)):
        assert bol_list[n], "No values exists"


def test_read_qc_table():
    # GIVEN the file exist
    # WHEN reading the file
    df = read_qc_table(HSMETRICS_QC_CHECK)

    # THEN check if the file contains any values
    bol_list = df.any().tolist()
    for n in range(len(bol_list)):
        assert bol_list[n], "No values exists"


def test_get_qc_value_and_qc_criteria_simple():
    # GIVEN the file exist
    hs_metrics_path = "tests/test_data/qc_files/multiqc_picard_HsMetrics.json"
    fold80 = "tests/test_data/qc_files/fold_80.txt"

    # WHEN reading the file
    qc_list = get_qc_value(read_json(hs_metrics_path))
    qc_criteria = check_qc_criteria_simple(qc_list, fold80)

    # THEN check if the qc-list in not empty and the qc-criteria is bool
    assert len(qc_list) > 0, "qc_list is empty"
    assert type(qc_criteria) == bool, "qc_criteria is not Boolean"


def test_get_bait_and_sample_name(tumor_normal_config):
    # GIVEN the file exists
    # WHEN reading the file
    bed = get_bait_name(tumor_normal_config)
    sample_name = get_sample_name(tumor_normal_config)

    # THEN check if bed is string format and if sample name exists
    assert isinstance(bed, str), "bed is not in string format"
    assert sample_name[0], "sample name doesn't exist"


def test_check_qc_criteria():
    # GIVEN following df and bed
    df_qc = read_qc_table(HSMETRICS_QC_CHECK)
    bed = "gmcksolid_4.1_hg19_design.bed"

    # WHEN reading the function
    criteria_df = get_qc_criteria(df_qc, bed)

    # THEN check if df has two columns
    nr_of_columns = list(criteria_df.columns)
    assert len(nr_of_columns) == 2, "number of columns != 2"


def test_check_qc_criteria_output_csv_and_qc(tmp_path, tumor_normal_config):
    # GIVEN following an output_path, an hs_metrics file, and a config_json with a matching bed name
    test_new_dir = tmp_path / "check_qc_results"
    test_new_dir.mkdir()
    output_path = test_new_dir / "output.csv"
    new_config_json_file = Path(
        test_new_dir / "new_config_tumor_normal.json").as_posix()

    with open(tumor_normal_config, 'r') as f:
        new_config_json = json.load(f)
    new_config_json["panel"][
        "capture_kit"] = "dummy_path/to/capture_kit/gmcksolid_4.1_hg19_design.bed"
    write_json(new_config_json, new_config_json_file)

    hs_metrics = "tests/test_data/qc_files/multiqc_picard_HsMetrics.json"

    qc_criteria_df = get_qc_criteria(read_qc_table(HSMETRICS_QC_CHECK),
                                     get_bait_name(new_config_json_file))
    hs_metrics_df = read_json(hs_metrics)
    sample_names = get_sample_name(new_config_json_file)

    # WHEN calling the functions
    extract_criteria = check_qc_criteria(qc_criteria_df, hs_metrics_df,
                                         sample_names[0], sample_names[1])
    write_output(extract_criteria, output_path)
    qc_check = failed_qc(extract_criteria, sample_names[0], sample_names[1])

    # THEN check if the output df has 3 indexes, if csv-file exists and if qc check is string
    nr_of_indexes = list(extract_criteria.index)
    assert len(nr_of_indexes) == 3
    assert output_path.exists(), "File doesn't exists"
    assert isinstance(qc_check, str), "qc_check not string"
