# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from pathlib import Path
import glob
import tempfile
import os
import logging


from BALSAMIC.utils.exc import BalsamicError


from BALSAMIC.constants.paths import BALSAMIC_DIR
from BALSAMIC.constants.analysis import FastqName, SampleType
from BALSAMIC.utils.io import write_finish_file
from BALSAMIC.utils.rule import get_threads, get_result_dir
from BALSAMIC.constants.workflow_params import WORKFLOW_PARAMS
from BALSAMIC.models.analysis import BalsamicWorkflowConfig, BalsamicConfigModel


# Initialize BalsamicConfigModel
balsamic = BalsamicConfigModel.parse_obj(config)

shell.prefix("set -eo pipefail; ")

localrules: all

LOG = logging.getLogger(__name__)
logging.getLogger("filelock").setLevel("WARN")

# parse parameters as constants to workflows
params = BalsamicWorkflowConfig.parse_obj(WORKFLOW_PARAMS)

# Get analysis dir
analysis_dir_home = balsamic.analysis.analysis_dir
# Get case id/name
case_id = balsamic.analysis.case_id
# Get result dir
result_dir = balsamic.analysis.result

analysis_dir = str(analysis_dir_home.joinpath(case_id)) + "/"
fastq_dir = str(result_dir.joinpath("fastq")) + "/"
qc_dir = str(result_dir.joinpath("qc")) + "/"
bam_dir = str(result_dir.joinpath("bam")) + "/"
cnv_dir = str(result_dir.joinpath("cnv")) + "/"

reffasta = config["reference"]["reference_genome"]
refgene_flat = config["reference"]["refgene_flat"]
access_5kb_hg19 = config["reference"]["access_regions"]
target_bed = config["panel"]["capture_kit"]
benchmark_dir = config["analysis"]["benchmark"]
version = config["analysis"]["pon_version"]

singularity_image = balsamic.singularity['image']

sample_names = balsamic.get_all_sample_names()

tmp_dir = os.path.join(result_dir, "tmp", "" )
Path.mkdir(Path(tmp_dir), parents=True, exist_ok=True)


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
pon_finish = expand(analysis_dir + "/" + "analysis_PON_finish")

config["rules"] = [
    "snakemake_rules/quality_control/fastp.rule",
    "snakemake_rules/align/sentieon_alignment.rule",
]

for r in config["rules"]:
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
        Path(singularity_image, "varcall_cnvkit.sif").as_posix()
    benchmark:
        Path(benchmark_dir, "cnvkit.targets.tsv").as_posix()
    shell:
        """
cnvkit.py target {input.target_bait} --annotate {input.refgene_flat} --split -o {output.target_bed};
cnvkit.py antitarget {input.target_bait} -g {input.access_bed} -o {output.offtarget_bed};
        """

rule create_coverage:
    input:
        bam = lambda wildcards: balsamic.get_final_bam_name(bam_dir = bam_dir, sample_name = wildcards.sample_name),
        target_bed = cnv_dir + "target.bed",
        antitarget_bed = cnv_dir + "antitarget.bed"
    output:
        target_cnn = cnv_dir + "{sample}.targetcoverage.cnn",
        antitarget_cnn = cnv_dir + "{sample}.antitargetcoverage.cnn"
    singularity:
        Path(singularity_image, "varcall_cnvkit.sif").as_posix()
    benchmark:
        Path(benchmark_dir, "cnvkit_{sample}.coverage.tsv").as_posix()
    shell:
        """
cnvkit.py coverage {input.bam} {input.target_bed} -o {output.target_cnn};
cnvkit.py coverage {input.bam} {input.antitarget_bed} -o {output.antitarget_cnn};
        """

rule create_reference:
    input:
        cnn = expand(cnv_dir + "{sample}.{prefix}coverage.cnn", sample=config["samples"], prefix=["target", "antitarget"]),
        ref = reffasta
    output:
        ref_cnn = pon_reference
    singularity:
        Path(singularity_image, "varcall_cnvkit.sif").as_posix()
    benchmark:
        Path(benchmark_dir, "cnvkit.reference.tsv").as_posix()
    shell:
        """
cnvkit.py reference {input.cnn} --fasta {input.ref} -o {output.ref_cnn} ;
        """
