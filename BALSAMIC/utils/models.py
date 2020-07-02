'''
Contains constants and models for analysis or filtering
'''
import BALSAMIC
import sys

from pydantic import BaseModel, ValidationError, validator, Field
from pydantic.types import DirectoryPath, FilePath
from typing import Optional, List, Dict

from datetime import datetime
from pathlib import Path

from BALSAMIC.utils.constants import CONDA_ENV_PATH, CONDA_ENV_YAML, RULE_DIRECTORY


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
    AD=VCFAttributes(tag_value=5, filter_name="balsamic_low_tumor_ad", field="INFO"),
    DP=VCFAttributes(tag_value=100, filter_name="balsamic_low_tumor_dp", field="INFO"),
    MQ=VCFAttributes(tag_value=50, filter_name="balsamic_low_mq", field="INFO"),
    AF_max=VCFAttributes(tag_value=1, filter_name="balsamic_af_one", field="INFO"),
    AF_min=VCFAttributes(tag_value=0.02, filter_name="balsamic_low_af", field="INFO"),
    varcaller_name="VarDict",
    filter_type="general",
    analysis_type="tumor_only",
    description="General purpose filters used for filtering VarDict")


class QCModel(BaseModel):
    """Contains settings for quality control and pre-processing"""

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
    """Contains analysis variables
    REQUIRED FIELDS: 
    case_id
    analysis_type
    sequencing_type,
    analysis_dir
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
    BALSAMIC_version: str = BALSAMIC.__version__
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
                    ) + f'_BALSAMIC_{BALSAMIC.__version__}_graph.pdf'

    @validator("config_creation_date")
    def datetime_as_string(cls, value):
        return datetime.now().strftime("%Y-%m-%d %H:%M")


class SampleInstanceModel(BaseModel):
    """Holds attributes for samples used in analysis"""
    file_prefix: str
    sample_type: str = Field(alias="type")
    readpair_suffix: List[str] = ["1", "2"]


class BioinfoToolsModel(BaseModel):
    """Holds versions of current bioinformatic tools used in analysis"""
    tabix: str
    bcftools: str
    fastqc: str
    manta: str
    picard: str
    bwa: str
    strelka: str
    gatk: str
    samtools: str
    sambamba: str
    vardic: Optional[str]
    cutadapt: Optional[str]


class PanelModel(BaseModel):
    """Holds attributes of PANEL BED file if provided"""

    capture_kit: Optional[str]
    chrom: Optional[List[str]]

    @validator("capture_kit")
    def path_as_abspath_str(cls, value):
        return str(Path(value).resolve())


class BalsamicConfigModel(BaseModel):
    """Summarizes config models in preparation for export """
	
    QC: QCModel
    vcf: VCFModel
    analysis: AnalysisModel
    samples: Dict[str, SampleInstanceModel]
    reference: Dict[str, Path]
    conda_env_yaml: FilePath = str(CONDA_ENV_YAML)
    rule_directory: DirectoryPath = str(RULE_DIRECTORY)
    singularity: FilePath
    bioinfo_tools: BioinfoToolsModel
    panel: Optional[PanelModel]

    @validator("reference")
    def abspath_as_str(cls, value):
        for k, v in value.items():
            value[k] = str(Path(v).resolve())
        return value

    @validator("singularity")
    def transform_path_to_dict(cls, value):
        return {"image": str(Path(value).resolve())}
