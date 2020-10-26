from BALSAMIC.utils.qc_check import read_hs_metrics
from BALSAMIC.utils.qc_check import read_qc_table
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
