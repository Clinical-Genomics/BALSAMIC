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
from BALSAMIC.constants.analysis import FastqName, SampleType, SequencingType
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


# Run information
reffasta: str = config_model.reference["reference_genome"]
refgene_flat: str = config_model.reference["refgene_flat"]
access_5kb_hg19: str = config_model.reference["access_regions"]
target_bed: str = config_model.panel.capture_kit
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

panel_name = os.path.split(target_bed)[1].replace('.bed','')

coverage_references = expand(cnv_dir + "{sample}.{cov}coverage.cnn", sample=config["samples"], cov=['target','antitarget'])
baited_beds = expand(cnv_dir + "{cov}.bed", cov=['target','antitarget'])
pon_reference = expand(cnv_dir + panel_name + "_CNVkit_PON_reference_" + version + ".cnn")
pon_finish = expand(analysis_dir + "analysis_PON_finish")

sequence_type = config['analysis']["sequencing_type"]
rules_to_include = []
if sequence_type == SequencingType.TARGETED:
    rules_to_include.append("snakemake_rules/quality_control/fastp_tga.rule")
else:
    rules_to_include.append("snakemake_rules/quality_control/fastp_wgs.rule")

rules_to_include.append("snakemake_rules/align/sentieon_alignment.rule")


for r in rules_to_include:
    include: Path(BALSAMIC_DIR, r).as_posix()

rule all:
    input:
        ref_cnn = pon_reference
    output:
        pon_finish_file = pon_finish
    run:
        write_finish_file(file_path=output.pon_finish_file)

rule create_target:
    input:
        target_bait = target_bed,
        refgene_flat = refgene_flat,
        access_bed = access_5kb_hg19
    output:
        target_bed = cnv_dir + "target.bed",
        offtarget_bed = cnv_dir + "antitarget.bed"
    singularity:
        Path(singularity_image, "cnvkit.sif").as_posix()
    benchmark:
        Path(benchmark_dir, "cnvkit.targets.tsv").as_posix()
    shell:
        """
cnvkit.py target {input.target_bait} --annotate {input.refgene_flat} --split -o {output.target_bed};
cnvkit.py antitarget {input.target_bait} -g {input.access_bed} -o {output.offtarget_bed};
        """

rule create_coverage:
    input:
        bam = lambda wildcards: config_model.get_final_bam_name(bam_dir = bam_dir, sample_name = wildcards.sample),
        target_bed = cnv_dir + "target.bed",
        antitarget_bed = cnv_dir + "antitarget.bed"
    output:
        target_cnn = cnv_dir + "{sample}.targetcoverage.cnn",
        antitarget_cnn = cnv_dir + "{sample}.antitargetcoverage.cnn"
    singularity:
        Path(singularity_image, "cnvkit.sif").as_posix()
    benchmark:
        Path(benchmark_dir, "cnvkit_{sample}.coverage.tsv").as_posix()
    shell:
        """
cnvkit.py coverage {input.bam} {input.target_bed} -o {output.target_cnn};
cnvkit.py coverage {input.bam} {input.antitarget_bed} -o {output.antitarget_cnn};
        """

rule create_reference:
    input:
        cnn = expand(cnv_dir + "{sample}.{prefix}coverage.cnn", sample=config_model.get_all_sample_names(), prefix=["target", "antitarget"]),
        ref = reffasta
    output:
        ref_cnn = pon_reference
    singularity:
        Path(singularity_image, "cnvkit.sif").as_posix()
    benchmark:
        Path(benchmark_dir, "cnvkit.reference.tsv").as_posix()
    shell:
        """
cnvkit.py reference {input.cnn} --fasta {input.ref} -o {output.ref_cnn} ;
        """
