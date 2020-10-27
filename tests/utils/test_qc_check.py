from BALSAMIC.utils.qc_check import read_hs_metrics, read_qc_table
from BALSAMIC.utils.qc_check import get_bait_name, get_sample_name, get_qc_criteria
from BALSAMIC.utils.constants import HSMETRICS_QC_CHECK


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
    try:
        open("/BALSAMIC/utils/constants.py")
    except IOError:
        print("File doesn't exist")

    # WHEN reading the file
    df = read_qc_table(HSMETRICS_QC_CHECK)

    # THEN check if the file contains any values
    bol_list = df.any().tolist()
    for n in range(len(bol_list)):
        assert bol_list[n], "No values exists"


def test_get_bait_name():
    # GIVEN the file exists
    config_file = "tests/test_data/qc_files/config_paired.json"

    # WHEN reading the file
    bed = get_bait_name(config_file)

    # THEN check if bed is string format
    assert isinstance(bed, str), "bed is not in string format"


def test_get_sample_name():
    # GIVEN the file exists
    config_file = "tests/test_data/qc_files/config_paired.json"

    # WHEN reading the file
    sample_name = get_sample_name(config_file)

    # THEN check if file name exsists
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