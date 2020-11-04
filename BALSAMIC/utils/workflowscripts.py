import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def get_file_contents(input_file, prefix_name):
    """ Reads the 2-column tsv file and returns file contents with header names.

    Arguments:
    input_file: Path to the TSV file 
    """

    input_df = pd.read_csv(input_file, sep='\t', header=None)
    input_df.columns = ['id', 'AF']
    input_df['method'] = prefix_name
    return input_df


def get_densityplot(input_file1, input_file2, prefix_name1, prefix_name2,
                    output_file):
    """ Reads two input dataframes and outputs a densityplot to the pdf file.
 
    Arguments:
    input_file1: Path to the TSV file1 (e.g: AFs calculated without consensuscall step) 
    input_file2: Path to the TSV file2 (e.g: AFs calculated after consensuscall step)
    prefix_name1: Label legend w.r.t file1 (e.g: "umiextract") 
    prefix_name2: Label legend w.r.t file2 (e.g: "consensuscall")
    output_file: Path to output filename with '.pdf' extn (e.g: "outplot.pdf")
    """

    dataframe1 = get_file_contents(input_file1, prefix_name1)
    dataframe2 = get_file_contents(input_file2, prefix_name2)
    sns.kdeplot(dataframe1['AF'],
                color='r',
                shade=True,
                label=set(dataframe1['method']))
    sns.kdeplot(dataframe2['AF'],
                color='b',
                shade=True,
                label=set(dataframe2['method']))
    plt.xlabel('Allelic Frequency (AF)')
    plt.ylabel('Probability Density')
    plt.savefig(output_file)
    return output_file
