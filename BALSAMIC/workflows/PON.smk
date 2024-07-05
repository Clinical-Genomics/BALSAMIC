# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import glob
import logging
import os
import tempfile
from pathlib import Path
from typing import Dict, List

from BALSAMIC.constants.analysis import FastqName, Gender, PONWorkflow, SampleType, SequencingType
from BALSAMIC.constants.paths import BALSAMIC_DIR
from BALSAMIC.constants.workflow_params import WORKFLOW_PARAMS, SLEEP_BEFORE_START
from BALSAMIC.models.config import ConfigModel
from BALSAMIC.models.params import BalsamicWorkflowConfig
from BALSAMIC.utils.exc import BalsamicError
from BALSAMIC.utils.io import write_finish_file
from BALSAMIC.utils.rule import get_fastp_parameters, get_result_dir, get_threads


# Initialize ConfigModel
config_model = ConfigModel.model_validate(config)

shell.prefix("set -eo pipefail; ")

localrules: all

LOG = logging.getLogger(__name__)

# parse parameters as constants to workflows
params = BalsamicWorkflowConfig.model_validate(WORKFLOW_PARAMS)

# Get case id/name
case_id: str = config_model.analysis.case_id
# Get analysis dir
analysis_dir_home: str = config_model.analysis.analysis_dir
analysis_dir: str = Path(analysis_dir_home, "analysis", case_id).as_posix() + "/"
# Get result dir
result_dir: str = Path(config_model.analysis.result).as_posix() + "/"

# Create a temporary directory with trailing /
tmp_dir: str = Path(result_dir, "tmp").as_posix() + "/"
Path.mkdir(Path(tmp_dir), parents=True, exist_ok=True)

# Directories
benchmark_dir: str = config_model.analysis.benchmark + "/"
fastq_dir: str = Path(result_dir, "fastq").as_posix() + "/"
bam_dir: str = Path(result_dir, "bam", "").as_posix() + "/"
cnv_dir: str = Path(result_dir, "cnv", "").as_posix() + "/"
qc_dir: str = Path(result_dir, "qc", "").as_posix() + "/"

# PON setting
pon_workflow: PONWorkflow = config_model.analysis.pon_workflow

# Run information
version: str = config_model.analysis.pon_version
singularity_image: str = config_model.singularity['image']
sample_names: List[str] = config_model.get_all_sample_names()

# Fastp parameters
fastp_parameters: Dict = get_fastp_parameters(config_model)

sequence_type = config['analysis']["sequencing_type"]
rules_to_include = []
rules_to_include.append("snakemake_rules/misc/sleep.rule")
if sequence_type == SequencingType.TARGETED:
    rules_to_include.append("snakemake_rules/quality_control/fastp_tga.rule")
else:
    rules_to_include.append("snakemake_rules/quality_control/fastp_wgs.rule")

rules_to_include.append("snakemake_rules/align/sentieon_alignment.rule")

if pon_workflow == PONWorkflow.CNVKIT:
    reffasta: str = config_model.reference["reference_genome"]
    refgene_flat: str = config_model.reference["refgene_flat"]
    access_5kb_hg19: str = config_model.reference["access_regions"]
    target_bed: str = config_model.panel.capture_kit
    panel_name = os.path.split(target_bed)[1].replace('.bed','')

    pon_reference = expand(cnv_dir + panel_name + "_CNVkit_PON_reference_" + version + ".cnn")
    rules_to_include.append("snakemake_rules/pon/cnvkit_create_pon.rule")

if pon_workflow in [PONWorkflow.GENS_MALE, PONWorkflow.GENS_FEMALE]:
    gender = Gender.MALE if pon_workflow == PONWorkflow.GENS_MALE else Gender.FEMALE

    pon_reference = expand(cnv_dir + "gens_pon_100bp.{gender}.{version}.hdf5", gender=gender, version=version)
    rules_to_include.append("snakemake_rules/variant_calling/gatk_read_counts.rule")
    rules_to_include.append("snakemake_rules/pon/gens_create_pon.rule")

pon_finish = expand(analysis_dir + "analysis_PON_finish")

for r in rules_to_include:
    include: Path(BALSAMIC_DIR, r).as_posix()

rule all:
    input:
        ref_cnn = pon_reference
    output:
        pon_finish_file = pon_finish
    run:
        write_finish_file(file_path=output.pon_finish_file)
