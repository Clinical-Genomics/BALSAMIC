import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def get_file_contents(input_file, prefix_name):
    """ Reads the 2-column tsv file and returns file contents with header names  """
    file = pd.read_csv(input_file, sep='\t', header=None)
    file.columns = ['id', 'AF']
    file['method'] = prefix_name
    return(file)

def get_densityplot(input_file1, input_file2, prefix_name1, prefix_name2, output_pdf):
    """ Reads two input dataframes and outputs a densityplot to the pdf file """
    dataframe1 = get_file_contents(input_file1,prefix_name1)
    dataframe2 = get_file_contents(input_file2, prefix_name2)
    sns.kdeplot(dataframe1['AF'], color='r', shade=True, label=set(dataframe1['method']))
    sns.kdeplot(dataframe2['AF'], color='b', shade=True, label=set(dataframe2['method']))
    plt.xlabel('Allelic Frequency (AF)')
    plt.ylabel('Probability Density')
    return(plt.savefig(output_pdf))
