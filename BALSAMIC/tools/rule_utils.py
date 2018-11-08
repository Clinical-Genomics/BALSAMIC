import os
import yaml


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

    result_dir = [
        config["analysis"]["analysis_dir"], config["analysis"]["sample_id"],
        config["analysis"]["result"]
    ]
    result_dir = os.path.normpath(os.path.join(*result_dir))
    return result_dir


def get_conda_env(yaml_file, pkg):
    """
    input: BALSAMIC_env.yaml file from BALSAMIC's installation, and a package's name
    output: string of conda env where packge is in
    """

    yaml_in = yaml.load(open(yaml_file))

    for s, pkgs in yaml_in.items():
        if pkg in pkgs:
            conda_env = s
            break

    return conda_env


def get_picard_mrkdup(config):
    """
    input: sample config file output from BALSAMIC
    output: mrkdup or rmdup strings
    """ 

    picard_str = "mrkdup"

    if "picard_rmdup" in config["QC"]:
        if config["QC"]["picard_rmdup"].upper() == "TRUE":
            picard_str = "rmdup"

    return picard_str

