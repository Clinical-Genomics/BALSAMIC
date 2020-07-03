from pathlib import Path
from datetime import datetime

from pydantic import BaseModel, ValidationError, validator, Field
from pydantic.types import DirectoryPath, FilePath
from typing import Optional, List, Dict

from BALSAMIC.utils.constants import CONDA_ENV_PATH, CONDA_ENV_YAML, RULE_DIRECTORY
from BALSAMIC import __version__ as BALSAMIC_version


class VCFAttributes(BaseModel):
    """General purpose filter to manage various VCF attributes

    This class handles three parameters for the purpose filtering variants
    based on a tag_values, filter_name, and which field in VCF.

    E.g. AD=VCFAttributes(tag_value=5, filter_name="balsamic_low_tumor_ad", field="INFO")
    A value of 5 from INFO field and filter_name will be balsamic_low_tumor_ad

    Attributes:
      tag_value: float
      filter_name: str
      field: str
    """

    tag_value: float
    filter_name: str
    field: str


class VarCallerFilter(BaseModel):
    """General purpose for variant caller filters

    This class handles attributes and filter for variant callers

    Attributes:
      AD: VCFAttributes (required); minimum allelic depth
      AF_min: VCFAttributes (optional); minimum allelic fraction
      AF_max: VCFAttributes (optional); maximum allelic fraction
      MQ: VCFAttributes (optional); minimum mapping quality
      DP: VCFAttributes (optional); minimum read depth
      varcaller_name: str (required); variant caller name
      filter_type: str (required); filter name for variant caller
      analysis_type: str (required); analysis type e.g. tumor_normal or tumor_only
      description: str (required); comment section for description
    """

    AD: VCFAttributes
    AF_min: Optional[VCFAttributes]
    AF_max: Optional[VCFAttributes]
    MQ: Optional[VCFAttributes]
    DP: VCFAttributes
    varcaller_name: str
    filter_type: str
    analysis_type: str
    description: str


VARDICT = VarCallerFilter(
    AD=VCFAttributes(tag_value=5,
                     filter_name="balsamic_low_tumor_ad",
                     field="INFO"),
    DP=VCFAttributes(tag_value=100,
                     filter_name="balsamic_low_tumor_dp",
                     field="INFO"),
    MQ=VCFAttributes(tag_value=50, filter_name="balsamic_low_mq",
                     field="INFO"),
    AF_max=VCFAttributes(tag_value=1,
                         filter_name="balsamic_af_one",
                         field="INFO"),
    AF_min=VCFAttributes(tag_value=0.02,
                         filter_name="balsamic_low_af",
                         field="INFO"),
    varcaller_name="VarDict",
    filter_type="general",
    analysis_type="tumor_only",
    description="General purpose filters used for filtering VarDict")


class QCModel(BaseModel):
    """Contains settings for quality control and pre-processing
        Attributes:
            picard_rmdup : Field(bool); whether duplicate removal is to be applied in the workflow
            adapter : Field(str(AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT)); adapter sequence to trim
            quality_trim : Field(bool); whether quality trimming it to be performed in the workflow
            adapter_trim : Field(bool); whether adapter trimming is to be performed in the workflow
            umi_trim : Field(bool); whether UMI trimming is to be performed in the workflow
            min_seq_length : Field(str(int)); minimum sequence length cutoff for reads
            umi_trim_length : Field(str(int)); length of UMI to be trimmed from reads
    
    """
    picard_rmdup: bool = False
    adapter: str = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
    quality_trim: bool = True
    adapter_trim: bool = False
    umi_trim: bool = False
    min_seq_length: int = 25
    umi_trim_length: int = 5

    @validator("min_seq_length", "umi_trim_length")
    def coerce_int_as_str(cls, value):
        return str(value)

    class Config:
        validate_all = True


class VCFModel(BaseModel):
    """Contains VCF config"""

    tnsnv: Dict[str, str] = {"mutation": "somatic", "type": "SNV"}
    manta: Dict[str, str] = {"mutation": "somatic", "type": "SV"}
    pindel: Dict[str, str] = {"mutation": "somatic", "type": "SV"}
    cnvkit: Dict[str, str] = {"mutation": "somatic", "type": "CNV"}
    mutect: Dict[str, str] = {"mutation": "somatic", "type": "SNV"}
    vardict: Dict[str, str] = {"mutation": "somatic", "type": "SNV"}
    strelka: Dict[str, str] = {"mutation": "somatic", "type": "SNV"}
    tnscope: Dict[str, str] = {"mutation": "somatic", "type": "SNV"}
    vcfmerge: Dict[str, str] = {"mutation": "somatic", "type": "SNV"}
    dnascope: Dict[str, str] = {"mutation": "germline", "type": "SNV"}
    tnhaplotyper: Dict[str, str] = {"mutation": "somatic", "type": "SNV"}
    manta_germline: Dict[str, str] = {"mutation": "germline", "type": "SV"}
    haplotypecaller: Dict[str, str] = {"mutation": "germline", "type": "SNV"}
    strelka_germline: Dict[str, str] = {"mutation": "germline", "type": "SNV"}


class AnalysisModel(BaseModel):
    """Pydantic model containing workflow variables

        Attributes:
            case_id : Field(required); string case identifier
            analysis_type : Field(required); string literal [single, paired]
                single : if only tumor samples are provided
                paired : if both tumor and normal samples are provided
            sequencing_type : Field(required); string literal [targeted, wgs]
                targeted : if capture kit was used to enrich specific genomic regions
                wgs : if whole genome sequencing was performed
            analysis_dir : Field(required); existing path where to save files

            fastq_path : Field(optional); Path where fastq files will be stored
            script : Field(optional); Path where snakemake scripts will be stored
            log : Field(optional); Path where logs will be saved
            result : Field(optional); Path where BALSAMIC output will be stored
            benchmark : Field(optional); Path where benchmark report will be stored
            dag : Field(optional); Path where DAG graph of workflow will be stored

            BALSAMIC_version  : Field(optional); Current version of BALSAMIC
            config_creation_date  : Field(optional); Timestamp when config was created
    """

    case_id: str
    analysis_type: str
    sequencing_type: str
    analysis_dir: DirectoryPath
    fastq_path: Optional[DirectoryPath]
    script: Optional[DirectoryPath]
    log: Optional[DirectoryPath]
    result: Optional[DirectoryPath]
    benchmark: Optional[DirectoryPath]
    dag: Optional[FilePath]
    BALSAMIC_version: str = BALSAMIC_version
    config_creation_date: Optional[str]

    class Config:
        validate_all = True

    @validator("analysis_dir")
    def dirpath_always_abspath(cls, value) -> str:
        return Path(value).resolve().as_posix()

    @validator("log")
    def parse_analysis_to_log_path(cls, value, values, **kwargs) -> str:
        return Path(values.get("analysis_dir"), values.get("case_id"),
                    "logs").as_posix() + "/"

    @validator("fastq_path")
    def parse_analysis_to_fastq_path(cls, value, values, **kwargs) -> str:
        return Path(values.get("analysis_dir"), values.get("case_id"),
                    "analysis", "fastq").as_posix() + "/"

    @validator("script")
    def parse_analysis_to_script_path(cls, value, values, **kwargs) -> str:
        return Path(values.get("analysis_dir"), values.get("case_id"),
                    "scripts").as_posix() + "/"

    @validator("result")
    def parse_analysis_to_result_path(cls, value, values, **kwargs) -> str:
        return Path(values.get("analysis_dir"), values.get("case_id"),
                    "analysis").as_posix()

    @validator("benchmark")
    def parse_analysis_to_benchmark_path(cls, value, values, **kwargs) -> str:
        return Path(values.get("analysis_dir"), values.get("case_id"),
                    "benchmarks").as_posix() + "/"

    @validator("dag")
    def parse_analysis_to_dag_path(cls, value, values, **kwargs) -> str:
        return Path(values.get("analysis_dir"), values.get("case_id"),
                    values.get("case_id")).as_posix(
                    ) + f'_BALSAMIC_{BALSAMIC_version}_graph.pdf'

    @validator("config_creation_date")
    def datetime_as_string(cls, value):
        return datetime.now().strftime("%Y-%m-%d %H:%M")


class SampleInstanceModel(BaseModel):
    """Holds attributes for samples used in analysis
    
        Attributes:
            file_prefix : Field(str); basename of sample pair
            sample_type : Field(str; alias=type); type of sample [tumor, normal]
            readpair_suffix : Field(List); currently always set to [1, 2]
        """

    file_prefix: str
    sample_type: str = Field(alias="type")
    readpair_suffix: List[str] = ["1", "2"]


class BioinfoToolsModel(BaseModel):
    """Holds versions of current bioinformatic tools used in analysis"""
    tabix: Optional[str]
    bcftools: Optional[str]
    fastqc: Optional[str]
    manta: Optional[str]
    picard: Optional[str]
    bwa: Optional[str]
    strelka: Optional[str]
    gatk: Optional[str]
    samtools: Optional[str]
    sambamba: Optional[str]
    vardic: Optional[str]
    cutadapt: Optional[str]


class PanelModel(BaseModel):
    """Holds attributes of PANEL BED file if provided
        Attributes:
            capture_kit : Field(str(Path)); string representation of path to PANEL BED file
            chrom : Field(list(str)); list of chromosomes in PANEL BED
    """

    capture_kit: Optional[str]
    chrom: Optional[List[str]]

    @validator("capture_kit")
    def path_as_abspath_str(cls, value):
        return Path(value).resolve().as_posix()


class BalsamicConfigModel(BaseModel):
    """Summarizes config models in preparation for export 
    
        Attributes:
            QC : Field(QCmodel); variables relevant for fastq preprocessing and QC
            vcf : Field(VCFmodel); variables relevand for variant calling pipeline
            samples : Field(Dict); dictionary containing samples submitted for analysis
            reference : Field(Dict); dictionary containign paths to reference genome files
            panel : Field(PanelModel(optional)); variables relevant to PANEL BED if capture kit is used
            bioinfo_tools : Field(BioinfoToolsModel); dictionary of bioinformatics software and their versions used for the analysis
            singularity : Field(Path); path to singularity container of BALSAMIC

            conda_env_yaml : Field(Path(CONVA_ENV_YAML)); path where Balsamic configs can be found
            rule_directory : Field(Path(RULE_DIRECTORY)); path where snakemake rules can be found

    """

    QC: QCModel
    vcf: VCFModel
    analysis: AnalysisModel
    samples: Dict[str, SampleInstanceModel]
    reference: Dict[str, Path]
    singularity: FilePath
    bioinfo_tools: BioinfoToolsModel
    conda_env_yaml: FilePath = CONDA_ENV_YAML
    rule_directory: DirectoryPath = RULE_DIRECTORY
    panel: Optional[PanelModel]

    @validator("reference")
    def abspath_as_str(cls, value):
        for k, v in value.items():
            value[k] = Path(v).resolve().as_posix()
        return value

    @validator("singularity")
    def transform_path_to_dict(cls, value):
        return {"image": Path(value).resolve().as_posix()}
