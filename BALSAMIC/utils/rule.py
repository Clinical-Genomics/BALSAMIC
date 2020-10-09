import re
import yaml
from pathlib import Path
import snakemake
from BALSAMIC.utils.cli import get_file_extension
from BALSAMIC.utils.cli import find_file_index
from BALSAMIC.utils.constants import (MUTATION_TYPE, MUTATION_CLASS,
                                      WORKFLOW_SOLUTION, ANALYSIS_TYPES)
from BALSAMIC.utils.exc import WorkflowRunError


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


def get_variant_callers(config, mutation_type: str, mutation_class: str,
                        analysis_type: str, workflow_solution: str):
    """ Get list of variant callers for a given list of input

    Args:
        config: A validated dictionary of case_config
        mutation_type: A mutation type string, e.g. SNV
        mutation_class: A mutation class string, e.g. somatic
        analysis_type: A analysis type string, e.g. paired
        workflow_solution: A workflow type string, e.g. BALSAMIC

    Returns:
        A list variant caller names extracted from config

    Raises:
        WorkflowRunError if mutation_type, mutation_class, analysis_type, or workflow_solution do not have valid value
    """
    valid_variant_callers = set()
    if mutation_type not in MUTATION_TYPE:
        raise WorkflowRunError(
            f"{mutation_type} is not a valid mutation type.")

    if workflow_solution not in WORKFLOW_SOLUTION:
        raise WorkflowRunError(
            f"{workflow_solution} is not a valid workflow solution.")

    if analysis_type not in ANALYSIS_TYPES:
        raise WorkflowRunError(
            f"{analysis_type} is not a valid analysis type.")

    if mutation_class not in MUTATION_CLASS:
        raise WorkflowRunError(
            f"{mutation_class} is not a valid mutation class.")

    for variant_caller_name, variant_caller_params in config["vcf"].items():
        if mutation_type in variant_caller_params.get(
                "type") and mutation_class in variant_caller_params.get(
                    "mutation") and analysis_type in variant_caller_params.get(
                        "analysis_type"
                    ) and workflow_solution in variant_caller_params.get(
                        "workflow_solution"):
            valid_variant_callers.add(variant_caller_name)
    return list(valid_variant_callers)


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


def get_rule_output(rules, rule_name, output_file_wildcards):
    """get list of existing output files from a given workflow

    Args:
        rule_names: rule_name to query from rules object
        rules: snakemake rules object
    
    Returns:
        output_files: list of tuples (file, file_index, rule_name, tags, id, file_extension) for rules
    """
    output_files = list()
    # Extract housekeeper tags from rule's params value
    housekeeper = getattr(rules, rule_name).params.housekeeper_id

    # Get temp_output files
    temp_files = getattr(rules, rule_name).rule.temp_output

    # Get list of named output from rule. e.g. output.vcf
    output_file_names = list(getattr(rules, rule_name).output._names.keys())

    for output_name in output_file_names:
        output_file = getattr(rules, rule_name).output[output_name]
        for file_wildcard_list in snakemake.utils.listfiles(output_file):
            file_to_store = file_wildcard_list[0]
            # Do not store file if it is a temp() output
            if file_to_store in temp_files:
                continue

            file_extension = get_file_extension(file_to_store)
            file_to_store_index = find_file_index(file_to_store)

            base_tags = list(file_wildcard_list[1])
            base_tags.append(output_name)

            delivery_id = get_delivery_id(
                id_candidate=housekeeper["id"],
                file_to_store=file_to_store,
                tags=base_tags,
                output_file_wildcards=output_file_wildcards)

            # Return empty string if delivery_id is not resolved.
            # This can happen when wildcard from one rule tries to match with a file
            # from another rule. example: vep_somatic might pick up ngs_filter_vardict files
            pattern = re.compile(r"{([^}\.[!:]+)")
            if pattern.findall(delivery_id):
                continue

            # Create a composit tag from housekeeper tag and named output
            composit_tag = "-".join([housekeeper["tags"], output_name])
            file_tags = base_tags + [composit_tag]

            # replace all instsances of "_" with "-", since housekeeper doesn't like _
            file_tags = [t.replace("_", "-") for t in file_tags]

            output_files.append(
                (file_to_store, file_to_store_index, rule_name,
                 ",".join(file_tags), delivery_id, file_extension))

            if file_to_store_index:
                for file_index in file_to_store_index:
                    # Create a composit tag from housekeeper tag and named output
                    composit_tag = "-".join(
                        [housekeeper["tags"], output_name, "index"])
                    file_index_tags = base_tags + [composit_tag]

                    # replace all instsances of "_" with "-", since housekeeper doesn't like _
                    file_index_tags = [
                        t.replace("_", "-") for t in file_index_tags
                    ]
                    output_files.append(
                        (file_index, str(), rule_name,
                         ",".join(file_index_tags), delivery_id,
                         get_file_extension(file_index)))

    return output_files


def get_delivery_id(id_candidate: str, file_to_store: str, tags: list,
                    output_file_wildcards: dict):
    """resolve delivery id from file_to_store, tags, and output_file_wildcards
  
    This function will get a filename, a list of tags, and an id_candidate. id_candidate should be form of a fstring.

    Args:
        id_candidate: a fstring format string. e.g. "{case_name}"
        file_to_store: a filename to search a resolved id
        tags: a list of tags with a resolve id in it
        output_file_wildcards: a dictionary of wildcards. Keys are wildcard names, and values are list of wildcard values
    
    Returns:
        delivery_id: a resolved id string. If it can't be resolved, it'll return the id_candidate value
    """

    delivery_id = id_candidate
    for resolved_id in snakemake.io.expand(id_candidate,
                                           **output_file_wildcards):
        if resolved_id in file_to_store and resolved_id in tags:
            delivery_id = resolved_id
            break

    return delivery_id


def get_reference_output_files(reference_files_dict: dict,
                               file_type: str) -> list:
    """ Returns list of files matching a file_type from reference files

    Args:
        reference_files_dict: A validated dict model from reference
        file_type: a file type string, e.g. vcf, fasta
  
    Returns:
        ref_vcf_list: list of file_type files that are found in reference_files_dict 
    """
    ref_vcf_list = []
    for reference_key, reference_item in reference_files_dict.items():
        if reference_item['file_type'] == file_type:
            ref_vcf_list.append(reference_item['output_file'])
    return ref_vcf_list
