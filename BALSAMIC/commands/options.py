"""Balsamic command options."""
import click

from BALSAMIC import __version__ as balsamic_version
from BALSAMIC.constants.analysis import (
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
    CLUSTER_MAIL_TYPES,
    CLUSTER_PROFILES,
    QOS,
    QOS_OPTIONS,
    ClusterProfile,
)
from BALSAMIC.constants.constants import LOG_LEVELS, LogLevel
from BALSAMIC.constants.rules import DELIVERY_RULES
from BALSAMIC.constants.workflow_params import VCF_DICT
from BALSAMIC.utils.cli import validate_cache_version, validate_exome_option

OPTION_ADAPTER_TRIM = click.option(
    "--adapter-trim/--no-adapter-trim",
    default=True,
    show_default=True,
    is_flag=True,
    help="Trim adapters from reads in FASTQ file",
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

OPTION_BENCHMARK = click.option(
    "--benchmark",
    default=False,
    is_flag=True,
    help="Profile slurm jobs. Make sure you have slurm profiler enabled in your HPC.",
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

OPTION_CLUSTER_CONFIG = click.option(
    "--cluster-config",
    type=click.Path(),
    help="Cluster configuration JSON file path",
)

OPTION_CLUSTER_MAIL = click.option(
    "--mail-user",
    type=click.STRING,
    help="User email to receive notifications from the cluster",
)

OPTION_CLUSTER_MAIL_TYPE = click.option(
    "--mail-type",
    type=click.Choice(CLUSTER_MAIL_TYPES),
    help="The mail type triggering cluster emails",
)

OPTION_CLUSTER_PROFILE = click.option(
    "-p",
    "--profile",
    show_default=True,
    default=ClusterProfile.SLURM,
    type=click.Choice(CLUSTER_PROFILES),
    help="Cluster profile to submit jobs",
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

OPTION_DISABLE_VARIANT_CALLER = click.option(
    "--disable-variant-caller",
    help=f"Run workflow with selected variant caller(s) disable. Use comma to remove multiple variant callers. Valid "
    f"values are: {list(VCF_DICT.keys())}",
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
    callback=validate_exome_option
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
    type=click.Choice([Gender.FEMALE, Gender.MALE]),
    default=Gender.FEMALE,
    show_default=True,
    help="Sample associated gender",
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

OPTION_QUALITY_TRIM = click.option(
    "--quality-trim/--no-quality-trim",
    default=True,
    show_default=True,
    is_flag=True,
    help="Trim low quality reads in FASTQ file",
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

OPTION_UMI = click.option(
    "--umi/--no-umi",
    default=True,
    show_default=True,
    is_flag=True,
    help="UMI processing steps for samples with UMI tags. For WGS cases, UMI is always disabled.",
)

OPTION_UMI_TRIM_LENGTH = click.option(
    "--umi-trim-length",
    default=5,
    show_default=True,
    type=click.INT,
    help="Trim N bases from reads in FASTQ file",
)
