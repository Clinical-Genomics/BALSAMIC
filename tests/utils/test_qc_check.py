from BALSAMIC.utils.qc_check import read_hs_metrics, read_qc_table, check_qc_criteria, write_output
from BALSAMIC.utils.qc_check import get_bait_name, get_sample_name, get_qc_criteria, failed_qc
from BALSAMIC.utils.constants import HSMETRICS_QC_CHECK
import os


def test_read_hs_metrics():
    # GIVEN the file exist
    hs_metrics_path = "tests/test_data/qc_files/multiqc_picard_HsMetrics.json"

    # WHEN reading the file
    df = read_hs_metrics(hs_metrics_path)

    # THEN check if the file contains any values
    bol_list = df.any().tolist()
    for n in range(len(bol_list)):
        assert bol_list[n], "No values exists"


def test_read_qc_table():
    # GIVEN the file exist
    from BALSAMIC.utils.constants import HSMETRICS_QC_CHECK

    # WHEN reading the file
    df = read_qc_table(HSMETRICS_QC_CHECK)

    # THEN check if the file contains any values
    bol_list = df.any().tolist()
    for n in range(len(bol_list)):
        assert bol_list[n], "No values exists"


def test_get_bait_and_sample_name():
    # GIVEN the file exists
    config_file = "tests/test_data/qc_files/case_config_tumor_normal.json"

    # WHEN reading the file
    bed = get_bait_name(config_file)
    sample_name = get_sample_name(config_file)

    # THEN check if bed is string format and if sample name exists
    assert isinstance(bed, str), "bed is not in string format"
    assert sample_name[0], "sample name doesn't exist"


def test_get_qc_criteria():
    # GIVEN following df and bed
    df_qc = read_qc_table(HSMETRICS_QC_CHECK)
    bed = "gmcksolid_4.1_hg19_design.bed"

    # WHEN reading the function
    criteria_df = get_qc_criteria(df_qc, bed)

    # THEN check if df has two columns
    nr_of_columns = list(criteria_df.columns)
    assert len(nr_of_columns) == 2, "number of columns != 2"


def test_check_qc_criteria_output_csv_and_qc():
    # GIVEN following variables
    config_file = "tests/test_data/qc_files/case_config_tumor_normal.json"
    hs_metrics = "tests/test_data/qc_files/multiqc_picard_HsMetrics.json"
    output_path = "tests/test_data/qc_files/output.csv"
    qc_criteria_df = get_qc_criteria(read_qc_table(HSMETRICS_QC_CHECK),
                                     get_bait_name(config_file))
    hs_metrics_df = read_hs_metrics(hs_metrics)
    sample_names = get_sample_name(config_file)

    # WHEN calling the functions
    extract_criteria = check_qc_criteria(qc_criteria_df, hs_metrics_df,
                                         sample_names[0], sample_names[1])
    write_output(extract_criteria, output_path)
    qc_check = failed_qc(extract_criteria, sample_names[0], sample_names[1])

    # THEN check if the output df has 3 indexes, if csv-file exists and if qc check is string
    nr_of_indexes = list(extract_criteria.index)
    assert len(nr_of_indexes) == 3, "number of indexes != 3"
    assert os.path.isfile(output_path), "File doesn't exists"
    assert isinstance(qc_check, str), "qc_check not string"
