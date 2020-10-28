import pandas as pd
import numpy as np
import json
import os
from BALSAMIC.utils.constants import HSMETRICS_QC_CHECK
from BALSAMIC.utils.rule import get_sample_type


def read_hs_metrics(hs_metrics_file: str):
    """Reads the HS metrics (json-format) and returns it as a DataFrame

    Args:
        hs_metrics_file: A string path to the file

    Returns:
        metrics_df: DataFrame

    """
    with open(hs_metrics_file):
        metrics_df = pd.read_json(hs_metrics_file)
    return metrics_df


def read_qc_table(qc_table: dict):
    """Reads the QC-table (json-format) and returns it as a DataFrame

    Args:
        qc_table: Dictionary imported from constants

    Returns:
        qc_df: DataFrame

    """
    qc_df = pd.DataFrame.from_dict(qc_table)
    return qc_df


def get_bait_name(input_config: str):
    """Get the bait name from case config

    Args:
        input_config: Path to config

    Returns:
        bait: string

    """
    with open(input_config) as f:
        load_config = json.load(f)

        # Read the config file and return the bait name from the json file
        bait = os.path.basename(load_config["panel"]["capture_kit"])

    return bait


def get_sample_name(input_config: str):
    """ Get the sample names from the config file

    Args:
        input_config: Path to config

    Returns:
        tumor, normal: string

    """
    with open(input_config) as f:
        load_config = json.load(f)

        # get_sample_type returns a list, extracting the sample name with [0]
        normal = get_sample_type(load_config["samples"], "normal")[0]
        tumor = get_sample_type(load_config["samples"], "tumor")[0]

    return normal, tumor


def get_qc_criteria(input_df: pd.DataFrame, bait: str) -> pd.DataFrame:
    """ Creates a new DataFrame with the QC criteria for only the desired bait set

    Args:
        input_df: qc table as DataFrame
        bait: desired bait as string

    Returns:
        qc_df: DataFrame

    """
    # Copy the desired columns
    qc_df = input_df[[bait, "METRIC_CRITERIA"]].copy()

    # Changing the column with the bait name
    qc_df = qc_df.rename(columns={bait: bait + "_criteria"})

    return qc_df


def check_qc_criteria(input_qc_df: pd.DataFrame,
                      input_hsmetrics_df: pd.DataFrame, normal_sample: str,
                      tumor_sample: str) -> pd.DataFrame:
    """ This function can be divided in different parts:
        1) Merging intersected values for the df with the desired QC criteria and bait set, with the HS Metrics df
        2) Creating new columns with the QC-differences from the QC criteria
        3) Setting QC flags
        4) Extract the columns with the QC flag as a new DataFrame

    Args:
        input_qc_df: DataFrame
        input_hsmetrics_df: DataFrame
        normal_sample: String
        tumor_sample: String

    Returns:
        qc_check_df: DataFrame

    """

    # 1) Merge the two df by col (axis = 1) for those rows that are shared (intersected) by passing join='inner'
    merged_df = pd.concat([input_hsmetrics_df, input_qc_df],
                          axis=1,
                          join='inner')
    column_header = list(merged_df.columns)

    # 2) Adding new col with the calculated difference in the qc values
    merged_df['qc_diff_' + normal_sample] = merged_df[
        column_header[2]] - merged_df[column_header[0]]
    merged_df['qc_diff_' + tumor_sample] = merged_df[
        column_header[2]] - merged_df[column_header[1]]

    # 3) Desired conditions for normal and tumor sample to pass. Two different conditions are required
    # since the conditions are different for the samples and should not overwrite each other.
    conditions_normal = [(merged_df['qc_diff_' + normal_sample] <= 0) &
                         (merged_df['METRIC_CRITERIA'] == 'gt'),
                         (merged_df['qc_diff_' + normal_sample] >= 0) &
                         (merged_df['METRIC_CRITERIA'] == 'lt')]

    conditions_tumor = [(merged_df['qc_diff_' + tumor_sample] <= 0) &
                        (merged_df['METRIC_CRITERIA'] == 'gt'),
                        (merged_df['qc_diff_' + tumor_sample] >= 0) &
                        (merged_df['METRIC_CRITERIA'] == 'lt')]

    # If above conditions are "True", set them as "pass"
    set_qc = ['Pass', 'Pass']

    # Adding new column with qc flag.
    merged_df['qc_check_' + normal_sample] = np.select(conditions_normal,
                                                       set_qc,
                                                       default="Fail")
    merged_df['qc_check_' + tumor_sample] = np.select(conditions_tumor,
                                                      set_qc,
                                                      default="Fail")

    # 4) create a new df and copy the desired columns (separated by ',').
    qc_check_df = merged_df[[
        'qc_check_' + normal_sample, 'qc_diff_' + normal_sample,
        'qc_check_' + tumor_sample, 'qc_diff_' + tumor_sample
    ]].copy()

    return qc_check_df


def failed_qc(input_df: pd.DataFrame, normal_sample: str,
              tumor_sample: str) -> pd.DataFrame:
    """ Outputs if the QC failed

    Args:
        input_df: DataFrame with qc parameters and qc differences
        normal_sample: String
        tumor_sample: String

    Returns:
        String

    """

    # copy the columns with qc criteria
    copy_qc_df = input_df[[
        'qc_check_' + normal_sample, 'qc_check_' + tumor_sample
    ]].copy()

    # Creating df which set "True" for "Fail" values and "False" for "Pass" values
    qc_boolean = copy_qc_df.isin(["Fail"])

    # Create a Series to check whether any element is set as True by .any() and convert to list with .tolist()
    boolean_check = qc_boolean.any().tolist()

    # Loop trough the list to check for True booleans which indicates for failed qc criteria.
    for n in range(len(boolean_check)):
        if boolean_check[n]:
            qc = "QC failed"
            return qc



def write_output(input_df: pd.DataFrame, output_path: str) -> pd.DataFrame:
    """ Outputs the QC parameters as csv-file

    Args:
        input_df: DataFrame with qc parameters and qc differences
        output_path: String with the desired output path

    Returns:
        CSV-file

    """

    output_df = input_df.to_csv(output_path, sep='\t')

    return output_df


def get_qc_check(hs_metrics, output, config):
    """ Runs all above functions to provide the desired outputs

    Args:
        hs_metrics: Path to hs_metrics file
        output: Path for output csv-file
        config: Path to case config

    Returns:
        CSV-file and prints if the QC-failed

    """
    # Read the HS metrics and qc table and convert to df
    hs_metrics_df = read_hs_metrics(hs_metrics)
    qc_table_df = read_qc_table(HSMETRICS_QC_CHECK)

    # Extract the bait name and create a new df with the desired qc criteria
    bait_set = get_bait_name(config)
    sample_names = get_sample_name(config)
    qc_criteria_df = get_qc_criteria(qc_table_df, bait_set)

    # Create a df with qc-flag for each criteria for each sample
    extract_criteria = check_qc_criteria(qc_criteria_df, hs_metrics_df,
                                         sample_names[0], sample_names[1])

    # Check if qc failed
    failed_qc(extract_criteria, sample_names[0], sample_names[1])

    write_output(extract_criteria, output)
