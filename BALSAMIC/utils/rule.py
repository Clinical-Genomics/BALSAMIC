import os
import re
import yaml
from pathlib import Path
import snakemake


def get_chrom(panelfile):
    """
    input: a panel bedfile
    output: list of chromosomes in the bedfile
    """

    lines = [line.rstrip('\n') for line in open(panelfile, 'r')]
    chrom = list(set([s.split('\t')[0] for s in lines]))
    return chrom


def get_vcf(config, var_caller, sample):
    """
    input: BALSAMIC config file
    output: retrieve list of vcf files
    """

    vcf = []
    for v in var_caller:
        for s in sample:
            vcf.append(config["vcf"][v]["type"] + "." +
                       config["vcf"][v]["mutation"] + "." + s + "." + v)
    return vcf


def get_sample_type(sample, bio_type):
    """
    input: sample dictionary from BALSAMIC's config file
    output: list of sample type id
    """

    type_id = []
    for sample_id in sample:
        if sample[sample_id]["type"] == bio_type:
            type_id.append(sample_id)
    return type_id


def get_result_dir(config):
    """
    input: sample config file from BALSAMIC
    output: string of result directory path
    """

    return config['analysis']['result']


def get_conda_env(yaml_file, pkg):
    """
    Retrieve conda environment for package from a predefined yaml file

    input: balsamic_env 
    output: string of conda env where packge is in
    """

    with open(yaml_file, 'r') as file_in:
        yaml_in = yaml.safe_load(file_in)

    conda_env_found = None

    for conda_env, pkgs in yaml_in.items():
        if pkg in pkgs:
            conda_env_found = conda_env
            break

    if conda_env_found is not None:
        return conda_env_found
    else:
        raise KeyError(f'Installed package {pkg} was not found in {yaml_file}')


def get_picard_mrkdup(config):
    """
    input: sample config file output from BALSAMIC
    output: mrkdup or rmdup strings
    """

    picard_str = "mrkdup"

    if "picard_rmdup" in config["QC"]:
        if config["QC"]["picard_rmdup"] == True:
            picard_str = "rmdup"

    return picard_str


def get_script_path(script_name: str):
    """
    Retrieves script path where name is matching {{script_name}}.
    """

    p = Path(__file__).parents[1]
    script_path = str(Path(p, 'assets/scripts', script_name))

    return script_path


def get_threads(cluster_config, rule_name='__default__'):
    """
    To retrieve threads from cluster config or return default value of 8
    """

    return cluster_config[rule_name]['n'] if rule_name in cluster_config else 8


def get_rule_output(rules, rule_names=[]):
    """get list of existing output files from a given workflow

    Args:
        rule_names: a list of rule names in the workflow. If no rules are given,
          then it will get all rules from the workflow.
        rules: snakemake rules object
    
    Returns:
        output_files: list of tuples (file_name, rule_name, wildcard) for rule_names
    """
#    with open('/home/hassan.foroughi/repos/BALSAMIC/run_tests/TN_panel/get_rule_output.log', 'w') as f:
#        print(rules, file = f)
    output_files = list()
    output_files.append(('output_file', 'rulename', 'wildcard_value'))
    if not rule_names:
        rule_names = vars(rules).keys()
    for my_rule in rule_names:
        for my_file in getattr(rules, my_rule).output:
            for file_wildcard_list in snakemake.utils.listfiles(my_file):
                output_files.append((file_wildcard_list[0], my_rule, list(file_wildcard_list[1])))
    
    return output_files


def get_rule_output_raw(rules, rule_names=[], output_file_wildcards={}):
    """get list of all possible output files from a given workflow

    Args:
        rule_names: a list of rule names in the workflow. If no rules are given,
          then it will get all rules from the workflow.
        rules: snakemake rules object
        output_file_wildcards: a dictionary with wildcards as keys and values as list of wildcard values 
    
    Returns:
        output_files: list of tuples (file_name, rule_name, wildcard) for rule_names
    """
    
    output_files_raw = list()
    output_files_raw.append(('output_file', 'rulename', 'wildcard_name'))
    wildcard_sets = set()
    if not rule_names:
        rule_names = vars(rules).keys()
    for my_rule in vars(rules).keys():
        if my_rule in rule_names:
            for my_file in getattr(rules, my_rule).output:
                pattern = re.compile(r"{([^}\.[!:]+)")
                wildcard_subset = dict()
                if pattern.findall(my_file):

                    wildcard_sets.update(pattern.findall(my_file))
                    for w in pattern.findall(my_file):
                        wildcard_subset[w] = output_file_wildcards[w]

                    for my_file_expanded in snakemake.io.expand(my_file,**output_file_wildcards):
                        output_files_raw.append((my_file_expanded, my_rule, list(wildcard_subset.keys())))
                else:
                    output_files_raw.append((my_file, my_rule, [])) 
    return output_files_raw
