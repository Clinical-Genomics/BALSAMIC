"""Balsamic command options."""

import click
from BALSAMIC import __version__ as balsamic_version
from BALSAMIC.constants.analysis import (
    SubmitSnakemake,
    ANALYSIS_WORKFLOWS,
    PON_WORKFLOWS,
    RUN_MODES,
    AnalysisWorkflow,
    Gender,
    PONWorkflow,
    RunMode,
)
from BALSAMIC.constants.cache import GENOME_VERSIONS, CacheVersion, GenomeVersion
from BALSAMIC.constants.cluster import (
    Partition,
    QOS,
    QOS_OPTIONS,
)
from BALSAMIC.constants.constants import LOG_LEVELS, LogLevel
from BALSAMIC.constants.rules import DELIVERY_RULES
from BALSAMIC.constants.paths import WORKFLOW_PROFILE, CACHE_PROFILE
from BALSAMIC.utils.cli import (
    validate_cache_version,
    validate_exome_option,
    validate_umi_min_reads,
)

OPTION_ANALYSIS_DIR = click.option(
    "--analysis-dir",
    type=click.Path(exists=True, resolve_path=True),
    required=True,
    help="Path to store the analysis results",
)

OPTION_ANALYSIS_WORKFLOW = click.option(
    "-w",
    "--analysis-workflow",
    default=AnalysisWorkflow.BALSAMIC,
    show_default=True,
    type=click.Choice(ANALYSIS_WORKFLOWS),
    help="Balsamic analysis workflow to be executed",
)
OPTION_ARTEFACT_SNV_OBSERVATIONS = click.option(
    "--artefact-snv-observations",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    help="VCF path of somatic SNVs called in high coverage normal samples (used in all workflows)",
)
OPTION_BACKGROUND_VARIANTS = click.option(
    "-b",
    "--background-variants",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    help="Background set of valid variants for UMI",
)

OPTION_BALSAMIC_CACHE = click.option(
    "--balsamic-cache",
    type=click.Path(exists=True, resolve_path=True),
    required=True,
    help="Path to BALSAMIC cache",
)

OPTION_CACHE_VERSION = click.option(
    "--cache-version",
    show_default=True,
    default=balsamic_version,
    type=click.STRING,
    callback=validate_cache_version,
    help=f"Cache version to be used for init or analysis. Use '{CacheVersion.DEVELOP}' or 'X.X.X'.",
)

OPTION_CADD_ANNOTATIONS = click.option(
    "--cadd-annotations",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    help="Path of CADD annotations",
)

OPTION_CANCER_GERMLINE_SNV_OBSERVATIONS = click.option(
    "--cancer-germline-snv-observations",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    help="VCF path of cancer germline SNV normal observations (WGS analysis workflow)",
)

OPTION_CANCER_SOMATIC_SNV_OBSERVATIONS = click.option(
    "--cancer-somatic-snv-observations",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    help="VCF path of cancer SNV tumor observations (WGS analysis workflow)",
)

OPTION_CANCER_SOMATIC_SNV_PANEL_OBSERVATIONS = click.option(
    "--cancer-somatic-snv-panel-observations",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    help="VCF path of cancer SNV tumor observations from matching gene panel",
)

OPTION_CANCER_SOMATIC_SV_OBSERVATIONS = click.option(
    "--cancer-somatic-sv-observations",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    help="VCF path of cancer SV observations (WGS analysis workflow)",
)

OPTION_CASE_ID = click.option(
    "--case-id",
    required=True,
    help="Sample ID for reporting, naming the analysis jobs, and analysis path",
)

OPTION_CLINICAL_SNV_OBSERVATIONS = click.option(
    "--clinical-snv-observations",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    help="VCF path of clinical SNV observations (WGS analysis workflow)",
)

OPTION_CLINICAL_SV_OBSERVATIONS = click.option(
    "--clinical-sv-observations",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    help="VCF path of clinical SV observations (WGS analysis workflow)",
)

OPTION_CLUSTER_ACCOUNT = click.option(
    "--account",
    type=click.STRING,
    help="Cluster account to run jobs",
)

OPTION_RUN_INTERACTIVELY = click.option(
    "--run-interactively",
    is_flag=True,
    default=False,
    help="Run Snakemake job submission interactively instead of submitting the submitter to cluster.",
)

OPTION_SOFT_FILTER_NORMAL = click.option(
    "--soft-filter-normal",
    is_flag=True,
    default=False,
    help="Flag to disable hard-filtering on presence of variants in matched normal sample",
)

OPTION_WORKFLOW_PARTITION = click.option(
    "--workflow-partition",
    show_default=True,
    type=click.STRING,
    default=Partition.CORE,
    help="Cluster node partition to run Snakemake jobs",
)

OPTION_HEADJOB_PARTITION = click.option(
    "--headjob-partition",
    type=str,
    required=False,
    default=None,
    help="Cluster node partition to run Snakemake head-job",
)

OPTION_MAX_RUN_HOURS = click.option(
    "--max-run-hours",
    required=False,
    show_default=True,
    default=SubmitSnakemake.MAX_RUN_HOURS,
    type=click.INT,
    help="The maximum number of hours that the sbatch script for snakemake is allowed to run on the cluster.",
)

OPTION_WORKFLOW_PROFILE = click.option(
    "--workflow-profile",
    show_default=True,
    type=click.Path(exists=True, resolve_path=True),
    default=WORKFLOW_PROFILE,
    help="Directory containing snakemake workflow profile specifying rule resources",
)

OPTION_CACHE_PROFILE = click.option(
    "--cache-profile",
    show_default=True,
    type=click.Path(exists=True, resolve_path=True),
    default=CACHE_PROFILE,
    help="Directory containing snakemake cache profile specifying rule resources for cache workflow",
)

OPTION_CLUSTER_QOS = click.option(
    "--qos",
    show_default=True,
    default=QOS.LOW,
    type=click.Choice(QOS_OPTIONS),
    help="QOS for cluster jobs",
)

OPTION_COSMIC_KEY = click.option(
    "-c",
    "--cosmic-key",
    required=False,
    type=click.STRING,
    help="Cosmic DB authentication key",
)

OPTION_DRAGEN = click.option(
    "--dragen",
    is_flag=True,
    default=False,
    help="Enable dragen variant caller",
)

OPTION_EXOME = click.option(
    "--exome",
    is_flag=True,
    default=False,
    help="Assign exome parameters to TGA workflow",
    callback=validate_exome_option,
)

OPTION_FASTQ_PATH = click.option(
    "--fastq-path",
    type=click.Path(exists=True, resolve_path=True),
    required=True,
    help="Path to directory containing unconcatenated FASTQ files",
)

OPTION_FORCE_ALL = click.option(
    "--force-all",
    show_default=True,
    default=False,
    is_flag=True,
    help="Force execution. This is equivalent to Snakemake --forceall.",
)

OPTION_GENDER = click.option(
    "--gender",
    required=False,
    type=click.Choice([Gender.FEMALE, Gender.MALE, Gender.UNKNOWN]),
    default=Gender.UNKNOWN,
    show_default=True,
    help="Case associated Gender",
)

OPTION_GENOME_VERSION = click.option(
    "-g",
    "--genome-version",
    show_default=True,
    default=GenomeVersion.HG19,
    type=click.Choice(GENOME_VERSIONS),
    help="Type and build version of the reference genome",
)

OPTION_GENOME_INTERVAL = click.option(
    "--genome-interval",
    required=False,
    type=click.Path(exists=True, resolve_path=True),
    help="Genome 100 bp interval-file (created with gatk PreprocessIntervals), used for GENS pre-processing.",
)

OPTION_GENS_COV_PON = click.option(
    "--gens-coverage-pon",
    required=False,
    type=click.Path(exists=True, resolve_path=True),
    help="GENS PON file, either male or female (created with gatk CreateReadCountPanelOfNormals), used for GENS pre-processing.",
)

OPTION_GNOMAD_AF5 = click.option(
    "--gnomad-min-af5",
    required=False,
    type=click.Path(exists=True, resolve_path=True),
    help="Gnomad VCF filtered to keep >= 0.05 AF, used for GENS pre-processing.",
)

OPTION_LOG_LEVEL = click.option(
    "--log-level",
    default=LogLevel.INFO,
    type=click.Choice(LOG_LEVELS),
    help="Logging level in terms of urgency",
    show_default=True,
)

OPTION_NORMAL_SAMPLE_NAME = click.option(
    "--normal-sample-name",
    required=False,
    type=click.STRING,
    help="Normal sample name",
)

OPTION_OUT_DIR = click.option(
    "-o",
    "--out-dir",
    required=True,
    type=click.Path(exists=True),
    help="Output directory for singularity containers and reference files",
)

OPTION_PANEL_BED = click.option(
    "-p",
    "--panel-bed",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    help="Panel bed file of target regions",
)

OPTION_PON_CNN = click.option(
    "--pon-cnn",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    help="Panel of normal reference (.cnn) for CNVkit",
)

OPTION_PON_WORKFLOW = click.option(
    "--pon-workflow",
    type=click.Choice(PON_WORKFLOWS),
    default=PONWorkflow.CNVKIT,
    required=True,
    help="Specify which PON to create.",
)

OPTION_PON_VERSION = click.option(
    "-v",
    "--version",
    default="v1",
    type=click.STRING,
    help="Version of the PON file to be generated",
)

OPTION_PRINT_FILES = click.option(
    "-p",
    "--print-files",
    is_flag=True,
    default=False,
    show_default=True,
    help="Print list of analysis files. Otherwise only final count will be printed.",
)

OPTION_QUIET = click.option(
    "-q",
    "--quiet",
    default=False,
    is_flag=True,
    help="Instruct Snakemake to not output any progress or rule information",
)

OPTION_RULES_TO_DELIVER = click.option(
    "-r",
    "--rules-to-deliver",
    multiple=True,
    default=DELIVERY_RULES,
    show_default=False,
    type=click.Choice(DELIVERY_RULES),
    help="Specify the rules to deliver. The delivery mode selected via the --delivery-mode option.",
)

OPTION_RUN_ANALYSIS = click.option(
    "-r",
    "--run-analysis",
    show_default=True,
    default=False,
    is_flag=True,
    help="Flag to run the actual analysis",
)

OPTION_RUN_MODE = click.option(
    "--run-mode",
    show_default=True,
    default=RunMode.CLUSTER,
    type=click.Choice(RUN_MODES),
    help="Run mode to execute Balsamic workflows",
)

OPTION_SAMPLE_CONFIG = click.option(
    "-s",
    "--sample-config",
    required=True,
    type=click.Path(),
    help="Sample configuration file",
)

OPTION_SENTIEON_INSTALL_DIR = click.option(
    "--sentieon-install-dir",
    type=click.Path(exists=True, resolve_path=True),
    required=True,
    help="Path to Sentieon install directory",
)

OPTION_SENTIEON_LICENSE = click.option(
    "--sentieon-license",
    required=True,
    type=click.STRING,
    help="Sentieon license in format IP:Port",
)

OPTION_SHOW_ONLY_MISSING_FILES = click.option(
    "-m",
    "--show-only-missing",
    is_flag=True,
    default=False,
    show_default=True,
    help="Only show missing analysis files.",
)

OPTION_SNAKEFILE = click.option(
    "-S",
    "--snakefile",
    type=click.Path(),
    help="Custom Snakefile for internal testing",
)

OPTION_SNAKEMAKE_OPT = click.option(
    "--snakemake-opt",
    multiple=True,
    help="Options to be passed to Snakemake",
)

OPTION_SWEGEN_SNV = click.option(
    "--swegen-snv",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    help="VCF path of Swegen SNV frequency database",
)

OPTION_SWEGEN_SV = click.option(
    "--swegen-sv",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    help="VCF path of Swegen SV frequency database",
)

OPTION_TUMOR_SAMPLE_NAME = click.option(
    "--tumor-sample-name",
    required=True,
    type=click.STRING,
    help="Tumor sample name",
)

OPTION_UMI_MIN_READS = click.option(
    "--umi-min-reads",
    type=click.STRING,
    callback=validate_umi_min_reads,
    help="Minimum raw reads supporting each UMI group. Format: 'x,y,z'.",
)

OPTION_ALLOWLIST_SNVS = click.option(
    "--allowlist-snvs",
    type=click.Path(exists=True, resolve_path=True),
    help="Path to allowlist file for SNVs, see read-the-docs for format.",
)
