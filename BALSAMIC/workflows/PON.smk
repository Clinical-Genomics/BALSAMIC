# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from pathlib import Path
import glob
import tempfile
import os
import logging

from typing import List, Dict
from BALSAMIC.utils.exc import BalsamicError


from BALSAMIC.constants.paths import BALSAMIC_DIR
from BALSAMIC.constants.analysis import FastqName, SampleType, SequencingType, PONWorkflow, Gender
from BALSAMIC.utils.io import write_finish_file
from BALSAMIC.utils.rule import get_fastp_parameters, get_threads, get_result_dir
from BALSAMIC.constants.workflow_params import WORKFLOW_PARAMS
from BALSAMIC.models.analysis import BalsamicWorkflowConfig, ConfigModel


# Initialize ConfigModel
config_model = ConfigModel.parse_obj(config)

shell.prefix("set -eo pipefail; ")

localrules: all

LOG = logging.getLogger(__name__)

# parse parameters as constants to workflows
params = BalsamicWorkflowConfig.parse_obj(WORKFLOW_PARAMS)

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
benchmark_dir: str = config_model.analysis.benchmark
fastq_dir: str = Path(result_dir, "fastq").as_posix() + "/"
bam_dir: str = Path(result_dir, "bam", "").as_posix() + "/"
cnv_dir: str = Path(result_dir, "cnv", "").as_posix() + "/"
qc_dir: str = Path(result_dir, "qc", "").as_posix() + "/"

# PON setting
pon_creation_type = config_model.analysis.pon_creation_type

# Run information
version: str = config_model.analysis.pon_version
singularity_image: str = config_model.singularity['image']
sample_names: List[str] = config_model.get_all_sample_names()

# Fastp parameters
fastp_parameters: Dict = get_fastp_parameters(config_model)

# Find and set Sentieon binary and license server from env variables
try:
    config["SENTIEON_LICENSE"] = os.environ["SENTIEON_LICENSE"]
    config["SENTIEON_INSTALL_DIR"] = os.environ["SENTIEON_INSTALL_DIR"]

    if os.getenv("SENTIEON_EXEC") is not None:
        config["SENTIEON_EXEC"] = os.environ["SENTIEON_EXEC"]
    else:
        config["SENTIEON_EXEC"] = Path(os.environ["SENTIEON_INSTALL_DIR"], "bin", "sentieon").as_posix()

except KeyError as error:
    LOG.error("Set environment variables SENTIEON_LICENSE, SENTIEON_INSTALL_DIR, SENTIEON_EXEC "
              "to run SENTIEON variant callers")
    raise BalsamicError

if not Path(config["SENTIEON_EXEC"]).exists():
    LOG.error("Sentieon executable not found {}".format(Path(config["SENTIEON_EXEC"]).as_posix()))
    raise BalsamicError

sequence_type = config['analysis']["sequencing_type"]
rules_to_include = []
if sequence_type == SequencingType.TARGETED:
    rules_to_include.append("snakemake_rules/quality_control/fastp_tga.rule")
else:
    rules_to_include.append("snakemake_rules/quality_control/fastp_wgs.rule")

rules_to_include.append("snakemake_rules/align/sentieon_alignment.rule")

if pon_creation_type == PONWorkflow.CNVKIT:
    reffasta: str = config_model.reference["reference_genome"]
    refgene_flat: str = config_model.reference["refgene_flat"]
    access_5kb_hg19: str = config_model.reference["access_regions"]
    target_bed: str = config_model.panel.capture_kit
    panel_name = os.path.split(target_bed)[1].replace('.bed','')
    coverage_references = expand(cnv_dir + "{sample}.{cov}coverage.cnn", sample=config["samples"], cov=['target',
                                                                                                      'antitarget'])
    baited_beds = expand(cnv_dir + "{cov}.bed",cov=['target', 'antitarget'])
    pon_reference = expand(cnv_dir + panel_name + "_CNVkit_PON_reference_" + version + ".cnn")
    rules_to_include.append("snakemake_rules/pon/cnvkit_create_pon.rule")

if pon_creation_type in [PONWorkflow.GENS_MALE, PONWorkflow.GENS_FEMALE]:
    if pon_creation_type == PONWorkflow.GENS_MALE:
        gender = Gender.MALE
    else:
        gender = Gender.FEMALE
    pon_reference = expand(cnv_dir + "balsamic_pon_100bp.{gender}.{version}.hdf5", gender=gender, version=version)
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

