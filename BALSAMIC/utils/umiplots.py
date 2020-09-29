import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def get_densityplot(input_file1, input_file2, output_file):
    umi= pd.read_csv(input_file1, sep='\t', header=None)
    umi.columns = ['id', 'AF']
    umi['method'] = 'umiextract'
    con = pd.read_csv(input_file2, sep='\t', header=None)
    con.columns = ['id', 'AF']
    con['method'] = 'consensuscall'
    sns.kdeplot(umi['AF'], color='r', shade=True, Label=set(umi['method']))
    sns.kdeplot(con['AF'], color='b', shade=True, Label=set(con['method']))
    plt.xlabel('Allelic Frequency (AF)')
    plt.ylabel('Probability Density')
    return (plt.savefig(output_file))
