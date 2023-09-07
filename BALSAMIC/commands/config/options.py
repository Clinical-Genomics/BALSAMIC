"""Balsamic config case command specific options."""
import click

from BALSAMIC import __version__ as balsamic_version
from BALSAMIC.constants.analysis import Gender, AnalysisWorkflow, ANALYSIS_WORKFLOWS
from BALSAMIC.constants.cache import GenomeVersion, GENOME_VERSIONS

OPTION_CASE_ID = click.option(
    "--case-id",
    required=True,
    help="Sample ID for reporting, naming the analysis jobs, and analysis path",
)

OPTION_GENDER = click.option(
    "--gender",
    required=False,
    type=click.Choice([Gender.FEMALE, Gender.MALE]),
    default=Gender.FEMALE.value,
    show_default=True,
    help="Sample associated gender",
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

OPTION_QUALITY_TRIM = click.option(
    "--quality-trim/--no-quality-trim",
    default=True,
    show_default=True,
    is_flag=True,
    help="Trim low quality reads in FASTQ file",
)

OPTION_ADAPTER_TRIM = click.option(
    "--adapter-trim/--no-adapter-trim",
    default=True,
    show_default=True,
    is_flag=True,
    help="Trim adapters from reads in FASTQ file",
)

OPTION_FASTQ_PATH = click.option(
    "--fastq-path",
    type=click.Path(exists=True, resolve_path=True),
    required=True,
    help="Path to directory containing unconcatenated FASTQ files",
)

OPTION_PANEL_BED = click.option(
    "-p",
    "--panel-bed",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    help="Panel bed file of target regions",
)

OPTION_BACKGROUND_VARIANTS = click.option(
    "-b",
    "--background-variants",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    help="Background set of valid variants for UMI",
)

OPTION_PON_CNN = click.option(
    "--pon-cnn",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    help="Panel of normal reference (.cnn) for CNVkit",
)

OPTION_BALSAMIC_CACHE = click.option(
    "--balsamic-cache",
    type=click.Path(exists=True, resolve_path=True),
    required=True,
    help="Path to BALSAMIC cache",
)

OPTION_CONTAINER_VERSION = click.option(
    "--container-version",
    show_default=True,
    default=balsamic_version,
    type=click.Choice(["develop", "master", balsamic_version]),
    help="Container for BALSAMIC version to download",
)

OPTION_ANALYSIS_DIR = click.option(
    "--analysis-dir",
    type=click.Path(exists=True, resolve_path=True),
    required=True,
    help="Path to store the analysis results",
)

OPTION_TUMOR_SAMPLE_NAME = click.option(
    "--tumor-sample-name",
    required=True,
    type=click.STRING,
    help="Tumor sample name",
)

OPTION_NORMAL_SAMPLE_NAME = click.option(
    "--normal-sample-name",
    required=False,
    type=click.STRING,
    help="Normal sample name",
)

OPTION_CADD_ANNOTATIONS = click.option(
    "--cadd-annotations",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    help="Path of CADD annotations",
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

OPTION_SWEGEN_SNV = click.option(
    "--swegen-snv",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    help="VCF path of Swegen SNV frequency database (WGS analysis workflow)",
)

OPTION_SWEGEN_SV = click.option(
    "--swegen-sv",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    help="VCF path of Swegen SV frequency database (WGS analysis workflow)",
)

OPTION_GENOME_VERSION = click.option(
    "-g",
    "--genome-version",
    default=GenomeVersion.HG19.value,
    type=click.Choice(GENOME_VERSIONS),
    help="Reference files genome version",
)

OPTION_ANALYSIS_WORKFLOW = click.option(
    "-w",
    "--analysis-workflow",
    default=AnalysisWorkflow.BALSAMIC.value,
    show_default=True,
    type=click.Choice(ANALYSIS_WORKFLOWS),
    help="Balsamic analysis workflow to be executed",
)

OPTION_PON_VERSION = click.option(
    "-v",
    "--version",
    default="v1",
    type=str,
    help="Version of the PON file to be generated",
)
