#! /usr/bin/python

import pandas as pd
import re
import numpy as np
import os
hs_metrics_path = '/Users/keyvan.elhami/Downloads/multiqc_picard_HsMetrics.json'
qc_table =   '/Users/keyvan.elhami/Downloads/qc_table4.json'
normal_sample='neatlyfastraven'
tumor_sample='easilyusefulorca'


def read_hs_metrics(hs_metrics_file: str):

    '''Reads the HS metrics (json-format) and returns it as a DataFrame

    Args:
        hs_metrics_file: A string path to the file

    Returns:
        metrics_df: DataFrame

    '''
    with open(hs_metrics_file) as f:
        metrics_df = pd.read_json(hs_metrics_file)
    return metrics_df


def read_qc_table(qc_table_file: str):

    '''Reads the QC-table (json-format) and returns it as a DataFrame

    Args:
        qc_table_file: A string path to the file

    Returns:
        qc_df: DataFrame

    '''
    with open(qc_table_file) as f:
        qc_df = pd.read_json(qc_table_file)
    return qc_df


def get_bait_name(input_df: pd.DataFrame) -> pd.DataFrame:

    '''Extracts the bait name from HS metrics DataFrame

    Args:
        input_df: HS metrics as DataFrame

    Returns:
        bait_name: string

    '''
    #"str.split" splits the bait name into a list with two elements.
    #exapnd = True will create new col in df with each element in the list as a value
    bait_set=input_df.loc['BAIT_SET'].str.split('_hg19_design.bed', expand = True)

    #dropping the last empty column which has 1 (int) as col name
    dropped_col_df = bait_set.drop(columns=1)

    #get the bed name and return it as a string
    bait_name = list(dropped_col_df.iloc[1])
    return str(bait_name[0])


def get_qc_criteria(input_df: pd.DataFrame, bait: str) -> pd.DataFrame:

    ''' Creates a new DataFrame with the QC criteria for only the deired bait set

    Args:
        input_df: qc table as DataFrame
        bait: desired bait as string

    Returns:
        qc_df: DataFrame

    '''
    qc_df = pd.DataFrame(data=input_df[bait])
    return qc_df


def check_qc_criteria(input_qc_df: pd.DataFrame, input_hsmetrics_df: pd.DataFrame) -> pd.DataFrame:

    ''' This function can be devided in different parts:
        1) Merging intersected values for the df with the desired QC criteria and bait set, with the HS Metrics df
        2) Creating new columns with the QC-differences from the QC criteria
        3) Setting QC flags
        4) Extract the columns with the QC flag as a new DataFrame

    Args:
        input_qc_df: DataFrame
        input_hsmetrics_df: DataFrame

    Returns:
        qc_check_df: DataFrame

    '''

    #1) Merge the two df by col (axis = 1) for those rows that are shared (intersected) by passing join='inner'
    merged_df = pd.concat([input_hsmetrics_df, input_qc_df], axis = 1, join='inner')

    column_header = list(merged_df.columns)

    #2) Adding new col with the calculated difference in the qc values
    merged_df['qc_diff_' + normal_sample] = merged_df[column_header[2]] - merged_df[column_header[0]]
    merged_df['qc_diff_' + tumor_sample] = merged_df[column_header[2]] - merged_df[column_header[1]]

    #3) Desired conditions for normal and tumor sample to pass. Two different conditions are required
    #since the conditions are different for the samples and should not overwrite each other.
    conditions_normal = [
        (merged_df['qc_diff_' + normal_sample] <= 0) & (merged_df['METRIC_CRITERIA'] == 'gt'),
        (merged_df['qc_diff_' + normal_sample] >= 0) & (merged_df['METRIC_CRITERIA'] == 'lt')]

    conditions_tumor = [
        (merged_df['qc_diff_' + tumor_sample] <= 0) & (merged_df['METRIC_CRITERIA'] == 'gt'),
        (merged_df['qc_diff_' + tumor_sample] >= 0) & (merged_df['METRIC_CRITERIA'] == 'lt')]

    #If above conditions are "True", set them as "pass"
    set_qc=['Pass', 'Pass']

    #Adding new column with qc flag.
    merged_df['qc_check_' + normal_sample] = np.select(conditions_normal, set_qc, default ="Fail")
    merged_df['qc_check_' + tumor_sample] = np.select(conditions_tumor, set_qc, default ="Fail")

    #4) create a new df and copy the desired columns (separated by ',').
    qc_check_df = merged_df[['qc_check_' + normal_sample, 'qc_diff_' + normal_sample,
                             'qc_check_' + tumor_sample, 'qc_diff_' + tumor_sample]].copy()

    return qc_check_df

@click.command()
@click.option('--hs_metrics', type = click.Path(exists=True), required = True, help = 'path to HS metrics for desired case')
@click.option('--qc_table', type = click.Path(exists=True), required = True, help = 'path to qc table with criteria' )

#The HS metrics and qc table provided in the command line will execute the main function.
def main(hs_metrics, qc_table):

    #Read the HS metrics and qc table and convert to df
    hs_metrics_df = read_hs_metrics(hs_metrics)
    qc_table_df = read_qc_table(qc_table)

    #Extract the bait name and create a new df with the desired qc criteria
    bait_set = get_bait_name(hs_metrics_df)
    qc_criteria_df = get_qc_criteria(qc_table_df, bait_set)

    #Create a df with qc-flag for each criteria for each sample
    check_qc_criteria(qc_criteria_df, hs_metrics_df)


if __name__ == '__main__':
    main()
