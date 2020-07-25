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


def get_rule_output(rules, rule_names: list, output_file_wildcards={}):
    """get list of existing output files from a given workflow

    Args:
        output_file_wildcards: a dictionary with wildcards as keys and values as list of wildcard values 
        rule_names: a list of rule names in the workflow. If no rules are given,
          then it will get all rules from the workflow.
        rules: snakemake rules object
    
    Returns:
        output_files: list of tuples (file_name, rule_name, wildcard) for rules
    """
    from BALSAMIC.utils.cli import find_file_index,get_file_extension
    output_files = list()
    output_files.append(('path', 'path_index', 'step', 'tag', 'id', 'format'))
    if not rule_names:
        rule_names = vars(rules).keys()
    rule_names = list(set(rule_names) & set(vars(rules).keys()))
    for my_rule in rule_names:
        delivery_id = getattr(rules, my_rule).params.housekeeper_id
        #delivery_id = snakemake.io.expand(delivery_id,**output_file_wildcards)
        for my_file in getattr(rules, my_rule).output:
#            pattern = re.compile(r"{([^}\.[!:]+)")
#            wildcard_sets = set()
#            if pattern.findall(my_file):
#                wildcard_sets.update(pattern.findall(my_file))
#            print(wildcard_sets)
            for file_wildcard_list in snakemake.utils.listfiles(my_file):
#                print(file_wildcard_list)
                file_to_store = file_wildcard_list[0]
                tags = list(file_wildcard_list[1]) 
                tags.extend(delivery_id["tags"])
                file_path_index = find_file_index(file_to_store)
                file_format = get_file_extension(my_file)
                if len(file_path_index) > 1:
                    LOG.warning("More than one index found for %s" % file_path_index)
                    LOG.warning("Taking %s index file" % file_path_index[0])
                
                file_path_index = file_path_index[0] if file_path_index else ""

                output_files.append((file_to_store, file_path_index, my_rule, ",".join(tags), ",".join(delivery_id["scope"]), file_format))
    
    return output_files


def get_rule_output_raw(rules, output_file_wildcards={}):
    """get list of all possible output files from a given workflow

    Args:
        rules: snakemake rules object
        output_file_wildcards: a dictionary with wildcards as keys and values as list of wildcard values 
    
    Returns:
        output_files: list of tuples (file_name, rule_name, wildcard) for rules
    """
    
    output_files_raw = list()
    output_files_raw.append(('output_file', 'rulename', 'wildcard_name'))
    wildcard_sets = set()
    for my_rule in vars(rules).keys():
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
