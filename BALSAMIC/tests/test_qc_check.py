from BALSAMIC.utils.qc_check import read_hs_metrics
from BALSAMIC.utils.qc_check import read_qc_table

def test_read_hs_metrics():
    #GIVEN the file exist
    hs_metrics_path = "/Users/keyvan.elhami/Downloads/multiqc_picard_HsMetrics.json"

    #WHEN reading the file
    df = read_hs_metrics(hs_metrics_path)

    #THEN check if the file contains any values
    bol_list = df.any().tolist()
    for n in range(len(bol_list)):
        assert bol_list[n], "No values exists"

def test_read_qc_table ():
    #GIVEN the file exist
    qc_table_path = "/Users/keyvan.elhami/Downloads/qc_table4.json"

    #WHEN reading the file
    df = read_qc_table(qc_table_path)

    #THEN check if the file contains any values
    bol_list = df.any().tolist()
    for n in range(len(bol_list)):
        assert bol_list[n], "No values exists"


