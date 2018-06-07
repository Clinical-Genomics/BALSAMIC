import os
import yaml

def get_vcf_files(config, vcf_type):
    vcf_files = []
    for var_caller in config["vcf"]:
        if config["vcf"][var_caller]["type"] == vcf_type:
            vcf_files.append(config["vcf"][var_caller]["merged"])
    return vcf_files


def get_sample_type(sample, bio_type):
    type_id = []
    for sample_id in sample:
        if sample[sample_id]["type"] == bio_type:
            type_id.append(sample_id)
    return type_id


def get_result_dir(config):
    result_dir = [
        config["analysis"]["analysis_dir"], config["analysis"]["sample_id"],
        config["analysis"]["result"]
    ]
    result_dir = os.path.normpath(os.path.join(*result_dir))
    return result_dir


def get_conda_env(yaml_file, pkg):
    
    yaml_in = yaml.load(open(yaml_file))

    for s, pkgs in yaml_in.items():
        if pkg in pkgs:
            conda_env = s
            break

    return conda_env

def get_picard_mrkdup(config):
    
    picard_str = "mrkdup"
    
    if "picard_rmdup" in config["QC"]:
        if config["QC"]["picard_rmdup"]:
            picard_str = "rmdup"

    return picard_str
