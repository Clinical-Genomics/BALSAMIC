# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
import os
import logging
import tempfile
import glob

from pathlib import Path
from yapf.yapflib.yapf_api import FormatFile

from snakemake.exceptions import RuleException, WorkflowError

from PyPDF2 import PdfFileMerger

from BALSAMIC.utils.exc import BalsamicError

from BALSAMIC.utils.cli import (check_executable, generate_h5)
from BALSAMIC.utils.io import write_json, read_yaml

from BALSAMIC.utils.models import VarCallerFilter, BalsamicWorkflowConfig

from BALSAMIC.utils.workflowscripts import plot_analysis

from BALSAMIC.utils.rule import (get_fastq_dict, get_variant_callers, get_rule_output, get_result_dir, get_vcf, get_picard_mrkdup,
                                 get_sample_type, get_threads, get_script_path, get_sequencing_type, get_capture_kit,
                                 get_clinical_snv_observations, get_clinical_sv_observations,get_swegen_snv,
                                 get_swegen_sv, dump_toml)

from BALSAMIC.constants.common import (SENTIEON_DNASCOPE, SENTIEON_TNSCOPE, RULE_DIRECTORY, MUTATION_TYPE)
from BALSAMIC.constants.variant_filters import (COMMON_SETTINGS, VARDICT_SETTINGS, SENTIEON_VARCALL_SETTINGS,
                                                SVDB_FILTER_SETTINGS)
from BALSAMIC.constants.workflow_params import (WORKFLOW_PARAMS, VARCALL_PARAMS)
from BALSAMIC.constants.workflow_rules import SNAKEMAKE_RULES

# Define Snakemake Rule Functions
def get_fastqpatterns():
    fastqpatterns = []
    for sample in fastq_dict:
        for fastqpattern in fastq_dict[sample]["fastqpair_patterns"]:
            fastqpatterns.append(fastqpattern)
    return fastqpatterns
def get_samplename(wcs):
    for sample in fastq_dict:
        if f"{wcs.fastqpattern}" in fastq_dict[sample]["fastqpair_patterns"]:
            return sample

def get_sampletype(wcs):
    for sample in fastq_dict:
        if f"{wcs.fastqpattern}" in fastq_dict[sample]["fastqpair_patterns"]:
            return fastq_dict[sample]["sampletype"]
def get_mapping(wcs):
    fastqpatterns = []
    for fastqpattern in fastq_dict[f"{wcs.sample}"]["fastqpair_patterns"]:
        fastqpatterns.append(fastqpattern)
    return expand(bam_dir + "{sample}_align_sort_{fastqpattern}.bam", sample=f"{wcs.sample}", fastqpattern=fastqpatterns)
def get_final_bam(wcs):
    if f"{wcs.sample_type}" == "tumor":
        samplename = tumor_sample
    else:
        samplename = normal_sample
    if config["analysis"]["sequencing_type"] == 'wgs':
        return bam_dir + "{sample}.dedup.realign.bam".format(sample=samplename)
    else:
        return bam_dir + "{tumor}.sorted.{picardstr}.bam".format(sample=samplename, picardstr=picarddup)
def get_final_tumor_bam(wcs=None):
    if config["analysis"]["sequencing_type"] == 'wgs':
        return bam_dir + "{tumor}.dedup.realign.bam".format(tumor=tumor_sample)
    else:
        return bam_dir + "{tumor}.sorted.{picardstr}.bam".format(tumor=tumor_sample, picardstr=picarddup)
def get_final_normal_bam(wcs=None):
    if config["analysis"]["sequencing_type"] == 'wgs':
        return bam_dir + "{normal}.dedup.realign.bam".format(normal=normal_sample)
    else:
        return bam_dir + "{normal}.sorted.{picardstr}.bam".format(normal=normal_sample, picardstr=picarddup)

shell.executable("/bin/bash")
shell.prefix("set -eo pipefail; ")

LOG = logging.getLogger(__name__)
logging.getLogger("filelock").setLevel("WARN")

# Create a temporary directory with trailing /
tmp_dir = os.path.join(get_result_dir(config), "tmp", "" )
Path.mkdir(Path(tmp_dir), parents=True, exist_ok=True)

# Set case id/name
case_id = config["analysis"]["case_id"]

# Directories
fastq_dir =  config["analysis"]["fastq_path"]
analysis_dir = config["analysis"]["analysis_dir"] + "/" +case_id + "/"
benchmark_dir = config["analysis"]["benchmark"]
fastqinput_dir = config["analysis"]["fastq_path"]
fastq_dir = get_result_dir(config) + "/fastq/"
bam_dir = get_result_dir(config) + "/bam/"
cnv_dir = get_result_dir(config) + "/cnv/"
fastqc_dir = get_result_dir(config) + "/fastqc/"
result_dir = get_result_dir(config) + "/"
vcf_dir = get_result_dir(config) + "/vcf/"
vep_dir = get_result_dir(config) + "/vep/"
qc_dir = get_result_dir(config) + "/qc/"
delivery_dir = get_result_dir(config) + "/delivery/"
umi_dir = get_result_dir(config) + "/umi/"
umi_qc_dir = qc_dir + "umi_qc/"
singularity_image = config['singularity']['image']


research_annotations = []
clinical_annotations = []
clinical_snv_obs = ""
swegen_snv = ""
clinical_sv = ""
swegen_sv = ""


# Prepare fastq_dict
fastq_dict = {}
for sample in config["samples"]:
    sampletype = config["samples"][sample]["type"]
    if sampletype == "tumor":
        tumor_sample = sample
        fastq_dict[tumor_sample] = get_fastq_dict(tumor_sample, fastqinput_dir)
        fastq_dict[tumor_sample]["sampletype"] = "TUMOR"
    else:
        normal_sample = sample
        fastq_dict[normal_sample] = get_fastq_dict(normal_sample,fastqinput_dir)
        fastq_dict[normal_sample]["sampletype"] = "NORMAL"


# vcfanno annotations
research_annotations.append( {
    'annotation': [{
    'file': Path(config["reference"]["gnomad_variant"]).as_posix(),
    'fields': ["AF", "AF_popmax"],
    'ops': ["self", "self"],
    'names': ["GNOMADAF", "GNOMADAF_popmax"]
    }]
}
)

research_annotations.append( {
    'annotation': [{
    'file': Path(config["reference"]["clinvar"]).as_posix(),
    'fields': ["CLNACC", "CLNREVSTAT", "CLNSIG", "ORIGIN", "CLNVC", "CLNVCSO"],
    'ops': ["self", "self", "self", "self", "self", "self"],
    'names': ["CLNACC", "CLNREVSTAT", "CLNSIG", "ORIGIN", "CLNVC", "CLNVCSO"]
    }]
}
)

if "swegen_snv_frequency" in config["reference"]:
    research_annotations.append( {
        'annotation': [{
            'file': get_swegen_snv(config),
            'fields': ["AF", "AC_Hom", "AC_Het", "AC_Hemi"],
            'ops': ["self", "self", "self","self"],
            'names': ["SWEGENAF", "SWEGENAAC_Hom", "SWEGENAAC_Het", "SWEGENAAC_Hemi"]
        }]
    }
    )

if "clinical_snv_observations" in config["reference"]:
    clinical_annotations.append( {
        'annotation': [{
            'file': get_clinical_snv_observations(config),
            'fields': ["Frq", "Obs", "Hom"],
            'ops': ["self", "self", "self"],
            'names': ["Frq", "Obs", "Hom"]
        }]
    }
    )
    clinical_snv_obs = get_clinical_snv_observations(config)


if "clinical_sv_observations" in config["reference"]:
    clinical_sv = get_clinical_sv_observations(config)


if "swegen_sv_frequency" in config["reference"]:
    swegen_sv = get_swegen_sv(config)


# picarddup flag
picarddup = get_picard_mrkdup(config)

# Varcaller filter settings
COMMON_FILTERS = VarCallerFilter.parse_obj(COMMON_SETTINGS)
VARDICT = VarCallerFilter.parse_obj(VARDICT_SETTINGS)
SENTIEON_CALLER = VarCallerFilter.parse_obj(SENTIEON_VARCALL_SETTINGS)
SVDB_FILTERS = VarCallerFilter.parse_obj(SVDB_FILTER_SETTINGS)

# parse parameters as constants to workflows
params = BalsamicWorkflowConfig.parse_obj(WORKFLOW_PARAMS)

# Capture kit name
if config["analysis"]["sequencing_type"] != "wgs":
    capture_kit = os.path.split(config["panel"]["capture_kit"])[1]


# explicitly check if cluster_config dict has zero keys.
if len(cluster_config.keys()) == 0:
    cluster_config = config

# Find and set Sentieon binary and license server from env variables
try:
    config["SENTIEON_LICENSE"] = os.environ["SENTIEON_LICENSE"]
    config["SENTIEON_INSTALL_DIR"] = os.environ["SENTIEON_INSTALL_DIR"]

    if os.getenv("SENTIEON_EXEC") is not None:
        config["SENTIEON_EXEC"] = os.environ["SENTIEON_EXEC"]
    else:
        config["SENTIEON_EXEC"] = Path(os.environ["SENTIEON_INSTALL_DIR"], "bin", "sentieon").as_posix()

    config["SENTIEON_TNSCOPE"] = SENTIEON_TNSCOPE
    config["SENTIEON_DNASCOPE"] = SENTIEON_DNASCOPE

except KeyError as error:
    LOG.error("Set environment variables SENTIEON_LICENSE, SENTIEON_INSTALL_DIR, SENTIEON_EXEC "
              "to run SENTIEON variant callers")
    raise BalsamicError

if not Path(config["SENTIEON_EXEC"]).exists():
    LOG.error("Sentieon executable not found {}".format(Path(config["SENTIEON_EXEC"]).as_posix()))
    raise BalsamicError

if "hg38" in config["reference"]["reference_genome"]:
    config["reference"]["genome_version"] = "hg38"
elif "canfam3" in config["reference"]["reference_genome"]:
    config["reference"]["genome_version"] = "canfam3"
    LOG.error("The main BALSAMIC workflow is not compatible with the canfam3 genome version "
             "use '--analysis-workflow balsamic-qc' instead")
    raise BalsamicError
else:
    config["reference"]["genome_version"] = "hg19"

LOG.info('Genome version set to %s', config["reference"]["genome_version"])


# Add normal sample if analysis is paired
germline_call_samples = ["tumor"]
if config['analysis']['analysis_type'] == "paired":
    germline_call_samples.append("normal")

# Create list of chromosomes in panel for panel only variant calling to be used in rules
if config["analysis"]["sequencing_type"] != "wgs":
    chromlist = config["panel"]["chrom"]

background_variant_file = ""
if "background_variants" in config:
    background_variant_file = config["background_variants"]

# Set temporary dir environment variable
os.environ["SENTIEON_TMPDIR"] = result_dir
os.environ['TMPDIR'] = get_result_dir(config)

# CNV report input files
cnv_data_paths = []
if config["analysis"]["sequencing_type"] == "wgs" and config['analysis']['analysis_type'] == "paired":
    cnv_data_paths.append(vcf_dir + "CNV.somatic." + config["analysis"]["case_id"] + ".ascat.samplestatistics.txt")
    cnv_data_paths.extend(expand(
        vcf_dir + "CNV.somatic." + config["analysis"]["case_id"] + ".ascat." + "{output_suffix}" + ".png",
        output_suffix=["ascatprofile", "rawprofile", "ASPCF", "tumor", "germline", "sunrise"]
    ))

if config["analysis"]["sequencing_type"] == "wgs" and config['analysis']['analysis_type'] == "single":
    cnv_data_paths.extend(expand(
        vcf_dir + "CNV.somatic." + config["analysis"]["case_id"] + ".cnvpytor." + "{output_suffix}" + ".png",
        output_suffix=["circular", "scatter"]
    ))

# Extract variant callers for the workflow
germline_caller = []
somatic_caller = []
somatic_caller_cnv = []
somatic_caller_sv = []
for m in MUTATION_TYPE:
    germline_caller_balsamic = get_variant_callers(config=config,
                                            analysis_type=config['analysis']['analysis_type'],
                                            workflow_solution="BALSAMIC",
                                            mutation_type=m,
                                            sequencing_type=config["analysis"]["sequencing_type"],
                                            mutation_class="germline")

    germline_caller_sentieon = get_variant_callers(config=config,
                                           analysis_type=config['analysis']['analysis_type'],
                                           workflow_solution="Sentieon",
                                           mutation_type=m,
                                           sequencing_type=config["analysis"]["sequencing_type"],
                                           mutation_class="germline")

    germline_caller = germline_caller + germline_caller_balsamic + germline_caller_sentieon


    somatic_caller_balsamic = get_variant_callers(config=config,
                                           analysis_type=config['analysis']['analysis_type'],
                                           workflow_solution="BALSAMIC",
                                           mutation_type=m,
                                           sequencing_type=config["analysis"]["sequencing_type"],
                                           mutation_class="somatic")

    somatic_caller_sentieon = get_variant_callers(config=config,
                                             analysis_type=config['analysis']['analysis_type'],
                                             workflow_solution="Sentieon",
                                             mutation_type=m,
                                             sequencing_type=config["analysis"]["sequencing_type"],
                                             mutation_class="somatic")

    somatic_caller_sentieon_umi = get_variant_callers(config=config,
                                             analysis_type=config['analysis']['analysis_type'],
                                             workflow_solution="Sentieon_umi",
                                             mutation_type=m,
                                             sequencing_type=config["analysis"]["sequencing_type"],
                                             mutation_class="somatic")
    somatic_caller = somatic_caller + somatic_caller_sentieon_umi + somatic_caller_balsamic + somatic_caller_sentieon

somatic_caller_sv = get_variant_callers(config=config,
                                            analysis_type=config['analysis']['analysis_type'],
                                            workflow_solution="BALSAMIC",
                                            mutation_type="SV",
                                            sequencing_type=config["analysis"]["sequencing_type"],
                                            mutation_class="somatic")

somatic_caller_cnv = get_variant_callers(config=config,
                                            analysis_type=config['analysis']['analysis_type'],
                                            workflow_solution="BALSAMIC",
                                            mutation_type="CNV",
                                            sequencing_type=config["analysis"]["sequencing_type"],
                                            mutation_class="somatic")
somatic_caller_sv.remove("svdb")
svdb_callers_prio = somatic_caller_sv + somatic_caller_cnv

for var_caller in svdb_callers_prio:
    if var_caller in somatic_caller:
        somatic_caller.remove(var_caller)

# Collect only snv callers for calculating tmb
somatic_caller_tmb = []
for ws in ["BALSAMIC","Sentieon","Sentieon_umi"]:
    somatic_caller_snv = get_variant_callers(config=config,
                                           analysis_type=config['analysis']['analysis_type'],
                                           workflow_solution=ws,
                                           mutation_type="SNV",
                                           sequencing_type=config["analysis"]["sequencing_type"],
                                           mutation_class="somatic")
    somatic_caller_tmb +=  somatic_caller_snv


# Remove variant callers from list of callers
if "disable_variant_caller" in config:
    variant_callers_to_remove = config["disable_variant_caller"].split(",")
    for var_caller in variant_callers_to_remove:
        if var_caller in somatic_caller:
            somatic_caller.remove(var_caller)
        if var_caller in germline_caller:
            germline_caller.remove(var_caller)

rules_to_include = []
analysis_type = config['analysis']["analysis_type"]
sequence_type = config['analysis']["sequencing_type"]

for sub,value in SNAKEMAKE_RULES.items():
  if sub in ["common", analysis_type + "_" + sequence_type]:
    for module_name,module_rules in value.items():
      rules_to_include.extend(module_rules)

if config["analysis"]["analysis_workflow"] == "balsamic":
    rules_to_include = [rule for rule in rules_to_include if "umi" not in rule]
    somatic_caller = [var_caller for var_caller in somatic_caller if "umi" not in var_caller]
    somatic_caller_tmb = [var_caller for var_caller in somatic_caller_tmb if "umi" not in var_caller]

# Adding code for testing, removing merge_bam from wgs analysis
if sequence_type == "wgs":
    if "snakemake_rules/variant_calling/mergetype_tumor.rule" in rules_to_include:
        rules_to_include.remove("snakemake_rules/variant_calling/mergetype_tumor.rule")
    if "snakemake_rules/variant_calling/mergetype_normal.rule" in rules_to_include:
        rules_to_include.remove("snakemake_rules/variant_calling/mergetype_normal.rule")

# Added code for testing, removing somalier because it had errors
if "snakemake_rules/quality_control/somalier.rule" in rules_to_include:
    rules_to_include.remove("snakemake_rules/quality_control/somalier.rule")

LOG.info(f"The following rules will be included in the workflow: {rules_to_include}")
LOG.info(f"The following Germline variant callers will be included in the workflow: {germline_caller}")
LOG.info(f"The following somatic variant callers will be included in the workflow: {somatic_caller}")


for r in rules_to_include:
    include: Path(RULE_DIRECTORY, r).as_posix()

# Define common and analysis specific outputs
quality_control_results = [
    os.path.join(qc_dir,case_id + "_metrics_deliverables.yaml"),
    os.path.join(qc_dir, "multiqc_report.html"),
    os.path.join(qc_dir, "multiqc_data/multiqc_data.json")
]

# Analysis results
analysis_specific_results = []

# Germline SNVs/SVs
analysis_specific_results.extend(
    expand(vep_dir + "{vcf}.vcf.gz", vcf=get_vcf(config, germline_caller, germline_call_samples))
)

# Germline SNVs specifically for genotype
if config["analysis"]["analysis_type"]=="paired":
    analysis_specific_results.append(vep_dir + "SNV.genotype.normal.dnascope.vcf.gz")

# Raw VCFs
analysis_specific_results.extend(
    expand(vcf_dir + "{vcf}.research.vcf.gz", vcf=get_vcf(config, somatic_caller, [case_id]))
)

# Filtered and passed post annotation research VCFs
analysis_specific_results.extend(
    expand(vep_dir + "{vcf}.research.filtered.pass.vcf.gz", vcf=get_vcf(config, somatic_caller, [case_id]))
)

# Filtered and passed post annotation clinical VCFs
analysis_specific_results.extend(
    expand(vep_dir + "{vcf}.clinical.filtered.pass.vcf.gz", vcf=get_vcf(config, somatic_caller, [case_id]))
)


# TMB
analysis_specific_results.extend(
    expand(vep_dir + "{vcf}.balsamic_stat", vcf=get_vcf(config, somatic_caller_tmb, [case_id]))
)

# WGS specific files
if config["analysis"]["sequencing_type"] == "wgs":
    # CNV report
    analysis_specific_results.append(vcf_dir + "CNV.somatic." + case_id + ".report.pdf"),

# TGA specific files
if config["analysis"]["sequencing_type"] != "wgs":
    # CNVkit
    analysis_specific_results.append(cnv_dir + "tumor.merged.cns")
    analysis_specific_results.extend(expand(cnv_dir + "tumor.merged-{plot}", plot=["diagram.pdf", "scatter.pdf"]))
    analysis_specific_results.append(cnv_dir + case_id +".gene_metrics")
    # vcf2cytosure
    analysis_specific_results.extend(expand(
        vcf_dir + "CNV.somatic.{case_name}.{var_caller}.vcf2cytosure.cgh",
        case_name=case_id,
        var_caller=["cnvkit"]
    ))
    # VarDict
    analysis_specific_results.extend(
        expand(vep_dir + "{vcf}.research.filtered.pass.ranked.vcf.gz", vcf=get_vcf(config, ["vardict"], [case_id]))
    )
    # UMI
    if config["analysis"]["analysis_workflow"]=="balsamic-umi":
        analysis_specific_results.extend(expand(umi_qc_dir + "{sample}.umi.mean_family_depth",sample=config["samples"]))
        if background_variant_file:
            analysis_specific_results.extend(
                expand(umi_qc_dir + "{case_name}.{var_caller}.AFtable.txt", case_name=case_id, var_caller=["tnscope_umi"])
        )

if config["analysis"]["sequencing_type"] == "wgs" and config['analysis']['analysis_type'] == "paired":
    analysis_specific_results.extend(
        expand(vcf_dir + "{vcf}.copynumber.txt.gz", vcf=get_vcf(config, ["ascat"], [case_id]))
    )
    analysis_specific_results.extend(
        expand(vcf_dir + "{vcf}.cov.gz", vcf=get_vcf(config,["dellycnv"],[case_id]))
    )
    analysis_specific_results.extend(expand(
        vcf_dir + "SV.somatic.{case_name}.{sample_type}.tiddit_cov.bed",
        case_name=case_id,
        sample_type=["tumor", "normal"]
    ))
    analysis_specific_results.extend(expand(
        vcf_dir + "CNV.somatic.{case_name}.{sample_type}.vcf2cytosure.cgh",
        case_name=case_id,
        sample_type=["tumor","normal"]
    ))

if config['analysis']['sequencing_type'] == "wgs" and config['analysis']['analysis_type'] == 'single':
    analysis_specific_results.extend(expand(
        vcf_dir + "CNV.somatic.{case_name}.{sample_type}.vcf2cytosure.cgh",
        case_name=case_id,
        sample_type=["tumor"]
    ))
    analysis_specific_results.extend(expand(
        vcf_dir + "SV.somatic.{case_name}.tumor.tiddit_cov.bed",
        case_name=case_id,
    ))

if config['analysis']['analysis_type'] == "single":
    analysis_specific_results.extend(
        expand(vcf_dir + "{vcf}.cov.gz",vcf=get_vcf(config,["dellycnv"],[case_id]))
    )

# Dragen
if config["analysis"]["sequencing_type"] == "wgs" and config['analysis']['analysis_type'] == "single":
    if "dragen" in config:
        analysis_specific_results.extend([
            Path(result_dir, "dragen", "SNV.somatic." + case_id + ".dragen_tumor.bam").as_posix(),
            Path(result_dir, "dragen", "SNV.somatic." + case_id + ".dragen.vcf.gz").as_posix()
        ])

LOG.info(f"Following outputs will be delivered {analysis_specific_results}")

if 'benchmark_plots' in config:
    log_dir = config["analysis"]["log"]
    if not check_executable("sh5util"):
        LOG.warning("sh5util executable does not exist. Won't be able to plot analysis")
    else:
        # Make individual plot per job
        for log_file in Path(log_dir).glob("*.err"):
            log_file_list = log_file.name.split(".")
            job_name = ".".join(log_file_list[0:4])
            job_id = log_file_list[4].split("_")[1]
            h5_file = generate_h5(job_name, job_id, log_file.parent)
            benchmark_plot = Path(benchmark_dir, job_name + ".pdf")

            log_file_plot = plot_analysis(log_file, h5_file, benchmark_plot)
            logging.debug("Plot file for {} available at: {}".format(log_file.as_posix(), log_file_plot))

        # Merge plots into one based on rule name
        for my_rule in vars(rules).keys():
            my_rule_pdf = PdfFileMerger()
            my_rule_plots = list()
            for plots in Path(benchmark_dir).glob(f"BALSAMIC*.{my_rule}.*.pdf"):
                my_rule_pdf.append(plots.as_posix())
                my_rule_plots.append(plots)
            my_rule_pdf.write(Path(benchmark_dir, my_rule+".pdf").as_posix())
            my_rule_pdf.close()

            # Delete previous plots after merging
            for plots in my_rule_plots:
                plots.unlink()

if 'delivery' in config:
    wildcard_dict = {
        "sample": list(config["samples"].keys())+["tumor", "normal"],
        "case_name": case_id,
        "allow_missing": True
    }

    if config['analysis']["analysis_type"] in ["paired", "single"]:
        wildcard_dict.update({
            "var_type": ["CNV", "SNV", "SV"],
            "var_class": ["somatic", "germline"],
            "var_caller": somatic_caller + germline_caller,
            "bedchrom": config["panel"]["chrom"] if "panel" in config else [],
        })

    if 'rules_to_deliver' in config:
        rules_to_deliver = config['rules_to_deliver'].split(",")
    else:
        rules_to_deliver = ['multiqc']

    output_files_ready = [('path', 'path_index', 'step', 'tag', 'id', 'format')]

    for my_rule in set(rules_to_deliver):
        try:
            housekeeper_id = getattr(rules, my_rule).params.housekeeper_id
        except (ValueError, AttributeError, RuleException, WorkflowError) as e:
            LOG.warning("Cannot deliver step (rule) {}: {}".format(my_rule, e))
            continue

        LOG.info("Delivering step (rule) {} {}.".format(my_rule, housekeeper_id))
        files_to_deliver = get_rule_output(rules=rules, rule_name=my_rule, output_file_wildcards=wildcard_dict)
        LOG.debug("The following files added to delivery: {}".format(files_to_deliver))
        output_files_ready.extend(files_to_deliver)

    output_files_ready = [dict(zip(output_files_ready[0], value)) for value in output_files_ready[1:]]
    delivery_ready = os.path.join(get_result_dir(config), "delivery_report", case_id + "_delivery_ready.hk")
    write_json(output_files_ready, delivery_ready)
    FormatFile(delivery_ready)

rule all:
    input:
        quality_control_results + analysis_specific_results
    output:
        finish_file = os.path.join(get_result_dir(config), "analysis_finish")
    params:
        tmp_dir = tmp_dir,
        case_name = config["analysis"]["case_id"],
    message:
        "Finalizing analysis for {params.case_name}",
    run:
        import datetime
        import shutil

        from BALSAMIC.utils.qc_metrics import validate_qc_metrics

        # Perform validation of extracted QC metrics
        try:
            validate_qc_metrics(read_yaml(input[0]))
        except ValueError as val_exc:
            LOG.error(val_exc)
            raise BalsamicError

        # Delete a temporal directory tree
        try:
            shutil.rmtree(params.tmp_dir)
        except OSError as e:
            print ("Error: %s - %s." % (e.filename, e.strerror))

        # Finish timestamp file
        with open(str(output.finish_file), mode="w") as finish_file:
            finish_file.write("%s\n" % datetime.datetime.now())
