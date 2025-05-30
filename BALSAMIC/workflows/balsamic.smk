# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
import glob
import logging
import os
import re
import tempfile
from pathlib import Path
from typing import Dict, List

from BALSAMIC.constants.constants import FileType
from BALSAMIC.constants.analysis import (
    AnalysisWorkflow,
    FastqName,
    MutationType,
    SampleType,
    SequencingType,
    AnalysisType,
    BioinfoTools)
from BALSAMIC.constants.paths import BALSAMIC_DIR
from BALSAMIC.constants.rules import SNAKEMAKE_RULES
from BALSAMIC.constants.variant_filters import (
    BaseSNVFilters,
    SVDB_FILTER_SETTINGS,
    MANTA_FILTER_SETTINGS,
    WgsSNVFilters,
    TgaSNVFilters,
    TgaUmiSNVFilters,
    get_tag_and_filtername,
)
from BALSAMIC.constants.workflow_params import (
    WORKFLOW_PARAMS,
    SLEEP_BEFORE_START,
)
from BALSAMIC.models.config import ConfigModel
from BALSAMIC.models.params import BalsamicWorkflowConfig, StructuralVariantFilters
from BALSAMIC.utils.cli import check_executable, generate_h5
from BALSAMIC.utils.exc import BalsamicError
from BALSAMIC.utils.io import read_yaml, write_finish_file, write_json
from BALSAMIC.utils.rule import (
    dump_toml,
    get_cancer_germline_snv_observations,
    get_cancer_somatic_snv_observations,
    get_capture_kit,
    get_clinical_snv_observations,
    get_clinical_sv_observations,
    get_fastp_parameters,
    get_pon_cnn,
    get_result_dir,
    get_rule_output,
    get_script_path,
    get_sequencing_type,
    get_somatic_sv_observations,
    get_swegen_snv,
    get_swegen_sv,
    get_threads,
    get_variant_callers,
    get_vcf,
)
from BALSAMIC.utils.workflowscripts import plot_analysis
from pypdf import PdfWriter
from snakemake.exceptions import RuleException, WorkflowError
from yapf.yapflib.yapf_api import FormatFile

# Initialize ConfigModel
config_model = ConfigModel.model_validate(config)

shell.executable("/bin/bash")
shell.prefix("set -eo pipefail; ")

LOG = logging.getLogger(__name__)

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
input_fastq_dir: str = config_model.analysis.fastq_path + "/"
benchmark_dir: str = config_model.analysis.benchmark + "/"
fastq_dir: str = Path(result_dir, "fastq").as_posix() + "/"
bam_dir: str = Path(result_dir, "bam").as_posix() + "/"
cnv_dir: str = Path(result_dir, "cnv").as_posix() + "/"
fastqc_dir: str = Path(result_dir, "fastqc").as_posix() + "/"
vcf_dir: str = Path(result_dir, "vcf").as_posix() + "/"
vep_dir: str = Path(result_dir, "vep").as_posix() + "/"
qc_dir: str = Path(result_dir, "qc").as_posix() + "/"
delivery_dir: str = Path(result_dir, "delivery").as_posix() + "/"
umi_dir: str = Path(result_dir, "umi").as_posix() + "/"
umi_qc_dir: str = Path(qc_dir, "umi_qc").as_posix() + "/"

# Annotations
research_annotations = []
clinical_annotations = []
artefact_snv_obs = ""
clinical_snv_obs = ""
cancer_germline_snv_obs = ""
cancer_somatic_snv_obs = ""
swegen_snv = ""
clinical_sv = ""
somatic_sv = ""
swegen_sv = ""

if config["analysis"]["sequencing_type"] != "wgs":
    pon_cnn: str = get_pon_cnn(config)

# Run information
singularity_image: str = config_model.singularity["image"]
sample_names: List[str] = config_model.get_all_sample_names()
tumor_sample: str = config_model.get_sample_name_by_type(SampleType.TUMOR)

analysis_type = config_model.analysis.analysis_type
sequencing_type = config_model.analysis.sequencing_type

if analysis_type == AnalysisType.PAIRED:
    normal_sample: str = config_model.get_sample_name_by_type(SampleType.NORMAL)

# Sample status to sampleID namemap
if analysis_type == AnalysisType.PAIRED:
    status_to_sample_id = (
        "TUMOR" + "\\\\t" + tumor_sample + "\\\\n" + "NORMAL" + "\\\\t" + normal_sample
    )
else:
    status_to_sample_id = "TUMOR" + "\\\\t" + tumor_sample

soft_filter_normal = config_model.analysis.soft_filter_normal

if config_model.panel:
    SNV_FILTERS = TgaSNVFilters
else:
    SNV_FILTERS = WgsSNVFilters

if sequencing_type == "targeted":
    exome = config_model.panel.exome
    snv_quality_filters = SNV_FILTERS.get_filters(category="quality", analysis_type=analysis_type, exome=exome)
    snv_research_filters = SNV_FILTERS.get_filters(category="research", analysis_type=analysis_type, exome=exome)
    snv_clinical_filters = SNV_FILTERS.get_filters(category="clinical", analysis_type=analysis_type, exome=exome)
else:
    snv_quality_filters = SNV_FILTERS.get_filters(category="quality", analysis_type=analysis_type)
    snv_research_filters = SNV_FILTERS.get_filters(category="research", analysis_type=analysis_type)
    snv_clinical_filters = SNV_FILTERS.get_filters(category="clinical", analysis_type=analysis_type)

if config_model.analysis.analysis_workflow == AnalysisWorkflow.BALSAMIC_UMI:
    SNV_FILTERS = TgaUmiSNVFilters
    umi_snv_quality_filters = SNV_FILTERS.get_filters(category="quality", analysis_type=analysis_type)
    umi_snv_research_filters = SNV_FILTERS.get_filters(category="research", analysis_type=analysis_type)
    umi_snv_clinical_filters = SNV_FILTERS.get_filters(category="clinical", analysis_type=analysis_type)



SVDB_FILTERS = StructuralVariantFilters.model_validate(SVDB_FILTER_SETTINGS)
MANTA_FILTERS = StructuralVariantFilters.model_validate(MANTA_FILTER_SETTINGS)

# Fastp parameters
fastp_parameters: Dict = get_fastp_parameters(config_model)

# parse parameters as constants to workflows
params = BalsamicWorkflowConfig.model_validate(WORKFLOW_PARAMS)

# Custom filter to set the minimum number of reads required to support each UMI tag group
if config_model.custom_filters and config_model.custom_filters.umi_min_reads:
    params.umiconsensuscall.filter_minreads = config_model.custom_filters.umi_min_reads

# vcfanno annotations
research_annotations.append(
    {
        "annotation": [
            {
                "file": Path(config["reference"]["gnomad_variant"]).as_posix(),
                "fields": ["AF", "AF_popmax"],
                "ops": ["self", "self"],
                "names": ["GNOMADAF", "GNOMADAF_popmax"],
            }
        ]
    }
)

research_annotations.append(
    {
        "annotation": [
            {
                "file": Path(config["reference"]["clinvar"]).as_posix(),
                "fields": [
                    "CLNVID",
                    "CLNREVSTAT",
                    "CLNSIG",
                    "ORIGIN",
                    "ONC",
                    "ONCDN",
                    "ONCREVSTAT",
                    "ONCDISDB",
                    "ONCCONF",
                    "CLNVC",
                    "CLNVCSO",
                ],
                "ops": ["self", "self", "self", "self", "self", "self", "self", "self", "self", "self", "self"],
                "names": [
                    "CLNVID",
                    "CLNREVSTAT",
                    "CLNSIG",
                    "ORIGIN",
                    "ONC",
                    "ONCDN",
                    "ONCREVSTAT",
                    "ONCDISDB",
                    "ONCCONF",
                    "CLNVC",
                    "CLNVCSO",
                ],
            }
        ]
    }
)

research_annotations.append(
    {
        "annotation": [
            {
                "file": Path(config["reference"]["cadd_snv"]).as_posix(),
                "names": ["CADD"],
                "ops": ["mean"],
                "columns": [6],
            }
        ]
    }
)


if "swegen_snv_frequency" in config["reference"]:
    research_annotations.append(
        {
            "annotation": [
                {
                    "file": get_swegen_snv(config),
                    "fields": ["AF", "AC_Hom", "AC_Het", "AC_Hemi"],
                    "ops": ["self", "self", "self", "self"],
                    "names": [
                        "SWEGENAF",
                        "SWEGENAAC_Hom",
                        "SWEGENAAC_Het",
                        "SWEGENAAC_Hemi",
                    ],
                }
            ]
        }
    )

if "artefact_snv_observations" in config["reference"]:
    clinical_annotations.append(
        {
            "annotation": [
                {
                    "file": Path(config["reference"]["artefact_snv_observations"]).as_posix(),
                    "fields": ["Frq", "Obs", "Hom"],
                    "ops": ["self", "self", "self"],
                    "names": ["ArtefactFrq", "ArtefactObs", "ArtefactHom"],
                }
            ]
        }
    )
    artefact_snv_obs: str = Path(config["reference"]["artefact_snv_observations"]).as_posix()

if "clinical_snv_observations" in config["reference"]:
    clinical_annotations.append(
        {
            "annotation": [
                {
                    "file": get_clinical_snv_observations(config),
                    "fields": ["Frq", "Obs", "Hom"],
                    "ops": ["self", "self", "self"],
                    "names": ["Frq", "Obs", "Hom"],
                }
            ]
        }
    )
    clinical_snv_obs: str = get_clinical_snv_observations(config)

if "cancer_germline_snv_observations" in config["reference"]:
    clinical_annotations.append(
        {
            "annotation": [
                {
                    "file": get_cancer_germline_snv_observations(config),
                    "fields": ["Frq", "Obs", "Hom"],
                    "ops": ["self", "self", "self"],
                    "names": [
                        "Cancer_Germline_Frq",
                        "Cancer_Germline_Obs",
                        "Cancer_Germline_Hom",
                    ],
                }
            ]
        }
    )
    cancer_germline_snv_obs: str = get_cancer_germline_snv_observations(config)

if "cancer_somatic_snv_observations" in config["reference"]:
    clinical_annotations.append(
        {
            "annotation": [
                {
                    "file": get_cancer_somatic_snv_observations(config),
                    "fields": ["Frq", "Obs", "Hom"],
                    "ops": ["self", "self", "self"],
                    "names": [
                        "Cancer_Somatic_Frq",
                        "Cancer_Somatic_Obs",
                        "Cancer_Somatic_Hom",
                    ],
                }
            ]
        }
    )
    cancer_somatic_snv_obs: str = get_cancer_somatic_snv_observations(config)

if "clinical_sv_observations" in config["reference"]:
    clinical_sv: str = get_clinical_sv_observations(config)

if "cancer_somatic_sv_observations" in config["reference"]:
    somatic_sv: str = get_somatic_sv_observations(config)

if "swegen_sv_frequency" in config["reference"]:
    swegen_sv: str = get_swegen_sv(config)


# Capture kit name
if config["analysis"]["sequencing_type"] != "wgs":
    capture_kit = os.path.split(config["panel"]["capture_kit"])[1]

# explicitly check if cluster_config dict has zero keys.
if len(cluster_config.keys()) == 0:
    cluster_config = config

if "hg38" in config["reference"]["reference_genome"]:
    config["reference"]["genome_version"] = "hg38"
elif "canfam3" in config["reference"]["reference_genome"]:
    config["reference"]["genome_version"] = "canfam3"
    LOG.error(
        "The main BALSAMIC workflow is not compatible with the canfam3 genome version "
        "use '--analysis-workflow balsamic-qc' instead"
    )
    raise BalsamicError
else:
    config["reference"]["genome_version"] = "hg19"

LOG.info("Genome version set to %s", config["reference"]["genome_version"])


# Add normal sample if analysis is paired
germline_call_samples = ["tumor"]
if config["analysis"]["analysis_type"] == "paired":
    germline_call_samples.append("normal")

# Create list of chromosomes in panel for panel only variant calling to be used in rules
if config["analysis"]["sequencing_type"] != "wgs":
    chromlist = config["panel"]["chrom"]

background_variant_file = ""
if "background_variants" in config:
    background_variant_file = config["background_variants"]

# Set temporary dir environment variable
os.environ["SENTIEON_TMPDIR"] = result_dir
os.environ["TMPDIR"] = get_result_dir(config)

# CNV report input files
cnv_report_paths = []
if config["analysis"]["sequencing_type"] == "wgs":
    if config["analysis"]["analysis_type"] == "paired":
        cnv_report_paths.append(
            f"{vcf_dir}MSI.somatic.{config['analysis']['case_id']}.msisensorpro.msi.pdf"
        )
        cnv_report_paths.append(
            f"{vcf_dir}CNV.somatic.{config['analysis']['case_id']}.ascat.samplestatistics.txt.pdf"
        )
        cnv_report_paths.extend(
            expand(
                f"{vcf_dir}CNV.somatic.{config['analysis']['case_id']}.ascat.{{output_suffix}}.png.pdf",
                output_suffix=[
                    "ascatprofile",
                    "rawprofile",
                    "ASPCF",
                    "tumor",
                    "germline",
                    "sunrise",
                ],
            )
        )
    else:
        cnv_report_paths.extend(
            expand(
                f"{vcf_dir}CNV.somatic.{config['analysis']['case_id']}.cnvpytor.{{output_suffix}}.png.pdf",
                output_suffix=["circular", "scatter"],
            )
        )
else:
    if config["analysis"]["analysis_type"] == "paired":
        cnv_report_paths.append(
            f"{vcf_dir}MSI.somatic.{config['analysis']['case_id']}.msisensorpro.msi.pdf"
        )
    cnv_report_paths.extend(
        expand(f"{cnv_dir}tumor.merged-{{plot}}.pdf", plot=["diagram", "scatter"])
    )
    cnv_report_paths.append(
        f"{cnv_dir}CNV.somatic.{config['analysis']['case_id']}.purecn.purity.csv.pdf"
    )


# Collect all rules to be run

rules_to_include = []

for sub, value in SNAKEMAKE_RULES.items():
    if sub in ["common", analysis_type + "_" + sequencing_type]:
        for module_name, module_rules in value.items():
            rules_to_include.extend(module_rules)

if config["analysis"]["analysis_workflow"] == "balsamic":
    rules_to_include = [rule for rule in rules_to_include if "umi" not in rule]

# Add rule for DRAGEN
if "dragen" in config:
    rules_to_include.append("snakemake_rules/concatenation/concatenation.rule")

# Add rule for GENS
if "gnomad_min_af5" in config["reference"]:
    rules_to_include.append("snakemake_rules/variant_calling/gens_preprocessing.rule")
if "gnomad_min_af5" in config["reference"] and sequencing_type == SequencingType.WGS:
    rules_to_include.append("snakemake_rules/variant_calling/gatk_read_counts.rule")

LOG.info(f"The following rules will be included in the workflow: {rules_to_include}")

# If workflow includes UMI filtered results, add these results
wf_solutions = ["BALSAMIC", "Sentieon"]
if config["analysis"]["analysis_workflow"] == "balsamic-umi":
    wf_solutions.append("Sentieon_umi")

# Extract variant callers for the workflow
germline_caller_snv = []
germline_caller_sv = []
germline_caller_cnv = []
somatic_caller_snv = []
somatic_caller_cnv = []
somatic_caller_sv = []

# Collect list of variant callers to be run
for ws in wf_solutions:
    somatic_caller_snv += get_variant_callers(
        config=config,
        analysis_type=config["analysis"]["analysis_type"],
        workflow_solution=ws,
        mutation_type="SNV",
        sequencing_type=config["analysis"]["sequencing_type"],
        mutation_class="somatic",
    )
    somatic_caller_sv += get_variant_callers(
        config=config,
        analysis_type=config["analysis"]["analysis_type"],
        workflow_solution=ws,
        mutation_type="SV",
        sequencing_type=config["analysis"]["sequencing_type"],
        mutation_class="somatic",
    )
    somatic_caller_cnv += get_variant_callers(
        config=config,
        analysis_type=config["analysis"]["analysis_type"],
        workflow_solution=ws,
        mutation_type="CNV",
        sequencing_type=config["analysis"]["sequencing_type"],
        mutation_class="somatic",
    )
    germline_caller_snv += get_variant_callers(
        config=config,
        analysis_type=config["analysis"]["analysis_type"],
        workflow_solution=ws,
        mutation_type="SNV",
        sequencing_type=config["analysis"]["sequencing_type"],
        mutation_class="germline",
    )
    germline_caller_sv += get_variant_callers(
        config=config,
        analysis_type=config["analysis"]["analysis_type"],
        workflow_solution=ws,
        mutation_type="SV",
        sequencing_type=config["analysis"]["sequencing_type"],
        mutation_class="germline",
    )
    germline_caller_cnv += get_variant_callers(
        config=config,
        analysis_type=config["analysis"]["analysis_type"],
        workflow_solution=ws,
        mutation_type="CNV",
        sequencing_type=config["analysis"]["sequencing_type"],
        mutation_class="germline",
    )

LOG.info(
    f"The following Somatic SNV variant callers will be included in the workflow: {somatic_caller_snv}"
)
LOG.info(
    f"The following Somatic SV variant callers will be included in the workflow: {somatic_caller_sv}"
)
LOG.info(
    f"The following Somatic CNV variant callers will be included in the workflow: {somatic_caller_cnv}"
)
LOG.info(
    f"The following Germline SNV variant callers will be included in the workflow: {germline_caller_snv}"
)
LOG.info(
    f"The following Germline SV variant callers will be included in the workflow: {germline_caller_sv}"
)
LOG.info(
    f"The following Germline CNV variant callers will be included in the workflow: {germline_caller_cnv}"
)

somatic_caller_sv.remove("svdb")
svdb_callers_prio = somatic_caller_sv + somatic_caller_cnv

somatic_caller = somatic_caller_snv + ["svdb"]
final_somatic_snv_caller = somatic_caller_snv
if config["analysis"]["sequencing_type"] != "wgs":
    remove_caller_list = ["tnscope", "vardict"]
    for remove_caller in remove_caller_list:
        if remove_caller in somatic_caller:
            somatic_caller.remove(remove_caller)
            final_somatic_snv_caller.remove(remove_caller)

# Define common and analysis specific outputs
quality_control_results = [
    Path(qc_dir, case_id + "_metrics_deliverables.yaml").as_posix(),
    Path(qc_dir, "multiqc_report.html").as_posix(),
    Path(qc_dir, "multiqc_data/multiqc_data.json").as_posix(),
]

# Analysis results
analysis_specific_results = []

# Germline SNVs/SVs
germline_caller = germline_caller_cnv + germline_caller_sv + germline_caller_snv
analysis_specific_results.extend(
    expand(
        vep_dir + "{vcf}.vcf.gz",
        vcf=get_vcf(config, germline_caller, germline_call_samples),
    )
)

for r in rules_to_include:
    include: Path(BALSAMIC_DIR, r).as_posix()


# Germline SNVs specifically for genotype
if config["analysis"]["analysis_type"] == "paired":
    analysis_specific_results.append(vep_dir + "SNV.genotype.normal.dnascope.vcf.gz")


# Raw VCFs
analysis_specific_results.extend(
    expand(
        vcf_dir + "{vcf}.research.vcf.gz",
        vcf=get_vcf(config, somatic_caller, [case_id]),
    )
)

# Filtered and passed post annotation research VCFs
analysis_specific_results.extend(
    expand(
        vep_dir + "{vcf}.research.filtered.pass.vcf.gz",
        vcf=get_vcf(config, somatic_caller, [case_id]),
    )
)

# Scored clinical SNV VCFs
analysis_specific_results.extend(
    expand(
        vep_dir + "{vcf}.clinical.scored.vcf.gz",
        vcf=get_vcf(config, final_somatic_snv_caller, [case_id]),
    )
)

# Filtered and passed post annotation clinical VCFs
analysis_specific_results.extend(
    expand(
        vep_dir + "{vcf}.clinical.filtered.pass.vcf.gz",
        vcf=get_vcf(config, somatic_caller, [case_id]),
    )
)


# TMB
somatic_caller.remove("svdb")
analysis_specific_results.extend(
    expand(
        vep_dir + "{vcf}.balsamic_stat",
        vcf=get_vcf(config, somatic_caller, [case_id]),
    )
)

# CNV report
analysis_specific_results.append(cnv_dir + "CNV.somatic." + case_id + ".report.pdf"),

# TGA specific files
if config["analysis"]["sequencing_type"] != "wgs":
    # CNVkit
    analysis_specific_results.append(cnv_dir + "tumor.merged.cns")
    analysis_specific_results.extend(
        expand(cnv_dir + "tumor.merged-{plot}", plot=["diagram.pdf", "scatter.pdf"])
    )
    analysis_specific_results.append(cnv_dir + case_id + ".gene_metrics")
    # vcf2cytosure
    analysis_specific_results.extend(expand(
        vcf_dir + "CNV.somatic.{case_name}.{var_caller}.vcf2cytosure.cgh",
        case_name=case_id,
        var_caller=["cnvkit"]
    ))
    # UMI
    if config["analysis"]["analysis_workflow"] == "balsamic-umi":
        analysis_specific_results.extend(expand(umi_qc_dir + "tumor.{sample}.umi.mean_family_depth", sample=tumor_sample))
        if config['analysis']['analysis_type'] == "paired":
            analysis_specific_results.extend(expand(umi_qc_dir + "normal.{sample}.umi.mean_family_depth",sample=normal_sample))
        if background_variant_file:
            analysis_specific_results.extend(
                expand(
                    umi_qc_dir + "{case_name}.{var_caller}.AFtable.txt",
                    case_name=case_id,
                    var_caller=["tnscope_umi"],
                )
            )


if (
    config["analysis"]["sequencing_type"] == "wgs"
    and config["analysis"]["analysis_type"] == "paired"
):
    analysis_specific_results.extend(
        expand(
            vcf_dir + "{vcf}.copynumber.txt.gz",
            vcf=get_vcf(config, ["ascat"], [case_id]),
        )
    )
    analysis_specific_results.extend(
        expand(vcf_dir + "{vcf}.cov.gz", vcf=get_vcf(config, ["dellycnv"], [case_id]))
    )
    analysis_specific_results.extend(
        expand(
            vcf_dir + "SV.somatic.{case_name}.{sample_type}.tiddit_cov.bed",
            case_name=case_id,
            sample_type=["tumor", "normal"],
        )
    )
    analysis_specific_results.extend(
        expand(
            vcf_dir + "CNV.somatic.{case_name}.{sample_type}.vcf2cytosure.cgh",
            case_name=case_id,
            sample_type=["tumor", "normal"],
        )
    )

if (
    config["analysis"]["sequencing_type"] == "wgs"
    and config["analysis"]["analysis_type"] == "single"
):
    analysis_specific_results.extend(
        expand(
            vcf_dir + "CNV.somatic.{case_name}.{sample_type}.vcf2cytosure.cgh",
            case_name=case_id,
            sample_type=["tumor"],
        )
    )
    analysis_specific_results.extend(
        expand(
            vcf_dir + "SV.somatic.{case_name}.tumor.tiddit_cov.bed",
            case_name=case_id,
        )
    )

if config["analysis"]["analysis_type"] == "single":
    analysis_specific_results.extend(
        expand(vcf_dir + "{vcf}.cov.gz", vcf=get_vcf(config, ["dellycnv"], [case_id]))
    )

# GENS Outputs
if "gnomad_min_af5" in config["reference"]:
    analysis_specific_results.extend(
        expand(
            cnv_dir + "{sample}.{gens_input}.bed.gz",
            sample=sample_names,
            gens_input=["cov", "baf"],
        )
    )


# Dragen
if (
    config["analysis"]["sequencing_type"] == "wgs"
    and config["analysis"]["analysis_type"] == "single"
):
    if "dragen" in config:
        analysis_specific_results.extend(
            [
                Path(
                    result_dir, "dragen", "SNV.somatic." + case_id + ".dragen_tumor.bam"
                ).as_posix(),
                Path(
                    result_dir, "dragen", "SNV.somatic." + case_id + ".dragen.vcf.gz"
                ).as_posix(),
            ]
        )

LOG.info(f"Following outputs will be delivered {analysis_specific_results}")

if "benchmark_plots" in config:
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
            logging.debug(
                "Plot file for {} available at: {}".format(
                    log_file.as_posix(), log_file_plot
                )
            )

        # Merge plots into one based on rule name
        for my_rule in vars(rules).keys():
            my_rule_pdf = PdfWriter()
            my_rule_plots = list()
            for plots in Path(benchmark_dir).glob(f"BALSAMIC*.{my_rule}.*.pdf"):
                my_rule_pdf.append(plots.as_posix())
                my_rule_plots.append(plots)
            my_rule_pdf.write(Path(benchmark_dir, my_rule + ".pdf").as_posix())
            my_rule_pdf.close()

            # Delete previous plots after merging
            for plots in my_rule_plots:
                plots.unlink()

if "delivery" in config:
    wildcard_dict = {
        "sample": sample_names,
        "case_name": case_id,
        "allow_missing": True,
    }

    if config["analysis"]["analysis_type"] in ["paired", "single"]:
        wildcard_dict.update(
            {
                "var_type": ["CNV", "SNV", "SV"],
                "var_class": ["somatic", "germline"],
                "var_caller": somatic_caller + germline_caller,
                "bedchrom": config["panel"]["chrom"] if "panel" in config else [],
            }
        )

    if "rules_to_deliver" in config:
        rules_to_deliver = config["rules_to_deliver"].split(",")
    else:
        rules_to_deliver = ["multiqc"]

    output_files_ready = [("path", "path_index", "step", "tag", "id", "format")]

    for my_rule in set(rules_to_deliver):
        try:
            housekeeper_id = getattr(rules, my_rule).params.housekeeper_id
        except (ValueError, AttributeError, RuleException, WorkflowError) as e:
            LOG.warning("Cannot deliver step (rule) {}: {}".format(my_rule, e))
            continue

        LOG.info("Delivering step (rule) {} {}.".format(my_rule, housekeeper_id))
        files_to_deliver = get_rule_output(
            rules=rules, rule_name=my_rule, output_file_wildcards=wildcard_dict
        )
        LOG.info("The following files added to delivery: {}".format(files_to_deliver))
        output_files_ready.extend(files_to_deliver)

    output_files_ready = [
        dict(zip(output_files_ready[0], value)) for value in output_files_ready[1:]
    ]
    delivery_ready = Path(
        get_result_dir(config), "delivery_report", case_id + "_delivery_ready.hk"
    ).as_posix()
    write_json(output_files_ready, delivery_ready)
    FormatFile(delivery_ready)


wildcard_constraints:
    sample="|".join(sample_names),


rule all:
    input:
        quality_control_results + analysis_specific_results,
    output:
        finish_file=Path(get_result_dir(config), "analysis_finish").as_posix(),
    params:
        tmp_dir=tmp_dir,
        case_name=config["analysis"]["case_id"],
    message:
        "Finalizing analysis for {params.case_name}"
    run:
        import datetime
        import shutil

        from BALSAMIC.utils.metrics import validate_qc_metrics

        # Perform validation of extracted QC metrics
        try:
            validate_qc_metrics(read_yaml(input[0]))
        except ValueError as val_exc:
            LOG.error(val_exc)
            raise BalsamicError

            # Remove temporary directory tree
        try:
            shutil.rmtree(params.tmp_dir)
        except OSError as e:
            print("Error: %s - %s." % (e.filename, e.strerror))

            # Finish timestamp file
        write_finish_file(file_path=output.finish_file)
