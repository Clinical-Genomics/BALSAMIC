import logging
import os
import re
from pathlib import Path
from typing import Dict

import snakemake
import toml
from BALSAMIC.constants.paths import SCRIPT_DIR

from BALSAMIC.constants.analysis import (
    AnalysisType,
    MutationOrigin,
    MutationType,
    SequencingType,
    WorkflowSolution,
)
from BALSAMIC.models.config import ConfigModel
from BALSAMIC.utils.cli import find_file_index, get_file_extension
from BALSAMIC.utils.exc import WorkflowRunError

LOG = logging.getLogger(__name__)


def get_vcf(config, var_caller, sample):
    """
    input: BALSAMIC config file
    output: retrieve list of vcf files
    """

    vcf = []
    for v in var_caller:
        for s in sample:
            vcf.append(
                config["vcf"][v]["mutation_type"]
                + "."
                + config["vcf"][v]["mutation"]
                + "."
                + s
                + "."
                + v
            )
    return vcf


def get_variant_callers(
    config,
    mutation_type: str,
    mutation_class: str,
    analysis_type: str,
    workflow_solution: str,
    sequencing_type: str,
):
    """Get list of variant callers for a given list of input

    Args:
        config: A validated dictionary of case_config
        mutation_type: A mutation type string, e.g. SNV
        mutation_class: A mutation class string, e.g. somatic
        analysis_type: A analysis type string, e.g. paired
        workflow_solution: A workflow type string, e.g. BALSAMIC
        sequencing_type: A sequencing type, e.g. wgs or targeted

    Returns:
        A list variant caller names extracted from config

    Raises:
        WorkflowRunError if values are not valid
    """

    valid_variant_callers = list()
    if mutation_type not in set(MutationType):
        raise WorkflowRunError(f"{mutation_type} is not a valid mutation type.")

    if workflow_solution not in set(WorkflowSolution):
        raise WorkflowRunError(f"{workflow_solution} is not a valid workflow solution.")

    if analysis_type not in set(AnalysisType):
        raise WorkflowRunError(f"{analysis_type} is not a valid analysis type.")

    if mutation_class not in set(MutationOrigin):
        raise WorkflowRunError(f"{mutation_class} is not a valid mutation class.")

    if sequencing_type not in set(SequencingType):
        raise WorkflowRunError(f"{sequencing_type} is not a valid sequencing type.")

    for variant_caller_name, variant_caller_params in config["vcf"].items():
        if (
            mutation_type in variant_caller_params.get("mutation_type")
            and mutation_class in variant_caller_params.get("mutation")
            and analysis_type in variant_caller_params.get("analysis_type")
            and workflow_solution in variant_caller_params.get("workflow_solution")
            and sequencing_type in variant_caller_params.get("sequencing_type")
        ):
            valid_variant_callers.append(variant_caller_name)
    return list(valid_variant_callers)


def get_sequencing_type(config):
    """
    input: sample config file from BALSAMIC
    output: sequencing type string ("targeted" or "wgs")
    """

    return config["analysis"]["sequencing_type"]


def get_analysis_type(config):
    """
    input: sample config file from BALSAMIC
    output: analysis type string ("paired" or "single")
    """

    return config["analysis"]["analysis_type"]


def get_capture_kit(config):
    """
    input: sample config file from BALSAMIC
    output: panel bed name
    """

    if config["analysis"]["sequencing_type"] != "wgs":
        return os.path.basename(config["panel"]["capture_kit"])
    else:
        return None


def get_sample_type_from_sample_name(config, sample_name):
    """
    input: case config file from BALSAMIC, and sample_name
    output: sample type
    """
    for sample in config["samples"]:
        if sample_name == sample["name"]:
            return sample["type"]


def get_result_dir(config):
    """
    input: sample config file from BALSAMIC
    output: string of result directory path
    """

    return config["analysis"]["result"]


def get_script_path(script_name: str) -> str:
    """Return the path to the script matching the file name."""
    return Path(SCRIPT_DIR, script_name).as_posix()


def get_threads(cluster_config, rule_name="__default__"):
    """
    To retrieve threads from cluster config or return default value of 8
    """

    return cluster_config[rule_name]["n"] if rule_name in cluster_config else 8


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
                LOG.debug(
                    "File is tagged as temporary file in the workflow: {}".format(
                        file_to_store
                    )
                )
                continue

            file_extension = get_file_extension(file_to_store)
            file_to_store_index = find_file_index(file_to_store)
            base_tags = list(file_wildcard_list[1])
            base_tags.append(output_name)

            delivery_id = get_delivery_id(
                id_candidate=housekeeper["id"],
                file_to_store=file_to_store,
                tags=base_tags,
                output_file_wildcards=output_file_wildcards,
            )

            # Return empty string if delivery_id is not resolved.
            # This can happen when wildcard from one rule tries to match with a file
            # from another rule. example: vep_somatic might pick up ngs_filter_vardict files
            pattern = re.compile(r"{([^}\.[!:]+)")
            if pattern.findall(delivery_id):
                LOG.error(
                    "Problem in pattern matching the following: {}".format(delivery_id)
                )
                continue

            # Create a composit tag from housekeeper tag and named output
            composit_tag = "-".join([housekeeper["tags"], output_name])
            file_tags = base_tags + [composit_tag]

            # replace all instances of "_" with "-", since housekeeper doesn't like _
            file_tags = [t.replace("_", "-") for t in file_tags]

            LOG.info("Found the following delivery id: {}".format(delivery_id))
            LOG.info("Found the following file to store: {}".format(file_to_store))
            LOG.info("Above file is in the following rule: {}".format(rule_name))
            output_files.append(
                (
                    file_to_store,
                    file_to_store_index,
                    rule_name,
                    ",".join(file_tags),
                    delivery_id,
                    file_extension,
                )
            )

            if file_to_store_index:
                for file_index in file_to_store_index:
                    # Create a composit tag from housekeeper tag and named output
                    composit_tag = "-".join([housekeeper["tags"], output_name, "index"])
                    file_index_tags = base_tags + [composit_tag]

                    # replace all instsances of "_" with "-", since housekeeper doesn't like _
                    file_index_tags = [t.replace("_", "-") for t in file_index_tags]
                    output_files.append(
                        (
                            file_index,
                            str(),
                            rule_name,
                            ",".join(file_index_tags),
                            delivery_id,
                            get_file_extension(file_index),
                        )
                    )

    return output_files


def get_delivery_id(
    id_candidate: str, file_to_store: str, tags: list, output_file_wildcards: dict
):
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
    for resolved_id in snakemake.io.expand(id_candidate, **output_file_wildcards):
        if resolved_id in file_to_store and resolved_id in tags:
            delivery_id = resolved_id
            break

    return delivery_id


def get_pon_cnn(config: dict) -> str:
    """Returns path of pon_cnn for TGA workflow

    Args:
        config: a config dictionary

    Returns:
       Returns path of pon_cnn generated by cnvkit

    """
    return config["panel"]["pon_cnn"] if "pon_cnn" in config["panel"] else " "


def get_clinical_snv_observations(config: dict) -> str:
    """Returns path for clinical snv observations

    Args:
        config: a config dictionary

    Returns:
        Path for clinical_snv_observations vcf file

    """
    return Path(config["reference"]["clinical_snv_observations"]).as_posix()


def get_cancer_germline_snv_observations(config: dict) -> str:
    """Returns path for cancer germline snv observations

    Args:
        config: a config dictionary

    Returns:
        Path for cancer-germline-snv-observations vcf file

    """
    return Path(config["reference"]["cancer_germline_snv_observations"]).as_posix()


def get_cancer_somatic_snv_observations(config: dict) -> str:
    """Returns path for cancer somatic snv observations

    Args:
        config: a config dictionary

    Returns:
        Path for cancer-somatic-snv-observations vcf file

    """
    return Path(config["reference"]["cancer_somatic_snv_observations"]).as_posix()


def get_swegen_snv(config: dict) -> str:
    """Returns path for swegen snv frequencies

    Args:
        config: a config dictionary

    Returns:
        Path for swegen_snv vcf file

    """
    return Path(config["reference"]["swegen_snv_frequency"]).as_posix()


def get_clinical_sv_observations(config: dict) -> str:
    """Returns path for clinical sv observations

    Args:
        config: a config dictionary

    Returns:
        Path for clinical_sv_observations vcf file

    """
    return Path(config["reference"]["clinical_sv_observations"]).as_posix()


def get_somatic_sv_observations(config: dict) -> str:
    """Returns path for somatic sv observations

    Args:
        config: a config dictionary

    Returns:
        Path for cancer_somatic_sv_observations vcf file

    """
    return Path(config["reference"]["cancer_somatic_sv_observations"]).as_posix()


def get_swegen_sv(config: dict) -> str:
    """Returns path for swegen sv frequencies

    Args:
        config: a config dictionary

    Returns:
        Path for swegen_sv vcf file

    """
    return Path(config["reference"]["swegen_sv_frequency"]).as_posix()


def dump_toml(annotations: list) -> str:
    """Returns list of converted annotation in toml format

    Args:
        annotations: a list of toml

    Returns:
        toml_annotation: a string of toml annotation resources
    """
    toml_annotations = ""
    for annotation in annotations:
        toml_annotations += toml.dumps(annotation)
    return toml_annotations


def get_fastp_parameters(config_model: ConfigModel) -> Dict:
    """Returns a dictionary with parameters for the fastp rules.

    Args:
        config_model: The case config json instantiated as the ConfigModel

    Returns:
        fastp_parameters_dict: Dictionary with 1 or 2 lists, depending on if sequencing type has UMIs or not.

    """
    fastp_parameters_dict = {}

    # Add quality trimming parameters
    fastp_trim_qual = [
        "--trim_tail1",
        "1",
        "--n_base_limit",
        config_model.QC.n_base_limit,
        "--length_required",
        config_model.QC.min_seq_length,
        "--low_complexity_filter",
        "--trim_poly_g",
        "--dont_eval_duplication",
    ]

    # Add adapter trimming parameters
    fastp_trim_adapter = [
        "--detect_adapter_for_pe",
        "--dont_eval_duplication",
    ]

    # Add unique parameters to TGA trimming which is done separately for adapter and quality trim
    if config_model.analysis.sequencing_type == SequencingType.TARGETED:
        fastp_trim_qual.append("--disable_adapter_trimming")
        fastp_trim_adapter.extend(
            [
                "--disable_quality_filtering",
                "--length_required",
                config_model.QC.min_seq_length,
            ]
        )

    fastp_parameters_dict["fastp_trim_qual"] = fastp_trim_qual
    fastp_parameters_dict["fastp_trim_adapter"] = fastp_trim_adapter

    return fastp_parameters_dict
