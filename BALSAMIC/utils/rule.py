import os
import yaml
from pathlib import Path


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
