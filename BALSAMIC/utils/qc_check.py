import pandas as pd
import re
import numpy as np
import os
import click


qc_list = ['MEAN_TARGET_COVERAGE', 'FOLD_80_BASE_PENALTY', 'PCT_OFF_BAIT']


def read_hs_metrics(hs_metrics_file):

    '''
    input: string
    output: DataFrame

    Reads the HS metrics (json-format) and returns it as a DataFrame
    '''

    with open(hs_metrics_file) as f:
        metrics_df = pd.read_json(hs_metrics_file)
    return metrics_df


def read_qc_table(qc_table_file):

    '''
    input: string
    output: DataFrame

    Reads the QC-table (json-format) and returns it as a DataFrame
    '''

    with open(qc_table_file) as f:
        qc_df = pd.read_json(qc_table_file)
    return qc_df


def get_bait_name(input_df):

    '''
    input: DataFrame
    output: string

    Extracts the bait name from HS metrics DataFrame
    '''

    #"str.split" splits the bait name into a list with two elements.
    #exapnd = True will create new col in df with each element in the list as a value
    bait_set=input_df.loc['BAIT_SET'].str.split('_hg19_design.bed', expand = True)

    #dropping the last empty column which has 1 (int) as col name
    dropped_col_df = bait_set.drop(columns=1)

    #get the bed name and return it as a string
    bait_name = list(dropped_col_df.iloc[1])
    return str(bait_name[0])


def get_qc_criteria(input_df, bait):

    '''
    input: DataFrame, string
    output: DataFrame

    Creates a new DataFrame with the QC criteria for only the deired bait set
    '''

    qc_df = pd.DataFrame(data=input_df[bait])
    return qc_df


def check_qc_criteria(input_qc_df, input_hsmetrics_df):

    '''
    input: DataFrame, DataFrame
    output: DataFrame

    This function can be devided in different parts:
    1) Merging intersected values for the df with the desired QC criteria and bait set, with the HS Metrics df
    2) Creating new columns with the QC-differences from the QC criteria
    3) Setting QC flags
    4) Extract the columns with the QC flag as a new DataFrame
    '''

    #merge the two df by col (axis = 1) for those rows that are shared (intersected) by passing join='inner'
    merged_df = pd.concat([input_qc_df, input_hsmetrics_df], axis = 1, join='inner')

    normal_sample='neatlyfastraven'
    tumor_sample='easilyusefulorca'
    column_header = list(merged_df.columns)

    #Adding new col with the calculated difference in the qc values
    merged_df['qc_diff_' + normal_sample] = merged_df[column_header[0]] - merged_df[column_header[1]]
    merged_df['qc_diff_' + tumor_sample] = merged_df[column_header[0]] - merged_df[column_header[2]]

    #Getting the mean cov value
    normal_meanCov = merged_df.loc['MEAN_TARGET_COVERAGE', 'qc_diff_' + normal_sample]
    tumor_meanCov = merged_df.loc['MEAN_TARGET_COVERAGE', 'qc_diff_' + tumor_sample]

    #If mean cov qc difference is negative, change it to a pos value and vice versa.
    #This is required since the remaining qc-values needs to be positive to pass.
    if normal_meanCov < 0:
        merged_df.set_value('MEAN_TARGET_COVERAGE', 'qc_diff_' + normal_sample, abs(normal_meanCov))
    else:
        merged_df.set_value('MEAN_TARGET_COVERAGE', 'qc_diff_' + normal_sample, normal_meanCov * -1)

    if tumor_meanCov < 0:
        merged_df.set_value('MEAN_TARGET_COVERAGE', 'qc_diff_' + tumor_sample, abs(tumor_meanCov))
    else:
        merged_df.set_value('MEAN_TARGET_COVERAGE', 'qc_diff_' + normal_sample, normal_meanCov * -1)

    #Setting qc flag Fail or Pass as a new column
    merged_df['qc_check_' + normal_sample] = np.where(merged_df['qc_diff_' + normal_sample] >= 0, 'Pass', 'Fail')
    merged_df['qc_check_' + tumor_sample] = np.where(merged_df['qc_diff_' + tumor_sample] >= 0, 'Pass', 'Fail')

    #create a new df and copy the desired columns (separated by ',').
    qc_check_df = merged_df[['qc_check_' + normal_sample, 'qc_check_' + tumor_sample]].copy()

    print (qc_check_df)
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
