import BALSAMIC
import sys
import pydantic

from pydantic import BaseModel, ValidationError, validator, Field
from pydantic.types import DirectoryPath, FilePath
from typing import Optional, List, Dict
from typing_extensions import Literal
from datetime import date, datetime, time, timedelta
from pathlib import Path

CONDA_ENV_PATH = Path(sys.modules["BALSAMIC"].__file__).parent.resolve() / "conda"
CONDA_ENV_YAML = Path(sys.modules["BALSAMIC"].__file__).parent.resolve() / "config" / "balsamic_env.yaml"
RULE_DIRECTORY = Path(sys.modules["BALSAMIC"].__file__).parent.resolve()


class QCModel(BaseModel):
    picard_rmdup : bool = False
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
    tnsnv: Dict[str, str] = {"mutation": "somatic", "type": "SNV"}
    manta : Dict[str, str] = {"mutation": "somatic", "type": "SV"}
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
    analysis_type: Literal["single", "paired"]
    sequencing_type : Literal["targeted", "wgs"]
    analysis_dir: Path
    fastq_path: Path = "analysis_dir/case_id/analysis/fastq_path"
    script: Path = "analysis_dir/case_id/scripts"
    log: Path = "analysis_dir/case_id/logs"
    result: Path = "analysis_dir/case_id/analysis"
    benchmark: Path = "analysis_dir/case_id/benchmarks"
    dag : Path = 'analysis_dir/case_id_BALSAMIC_v.v.v_graph.pdf'
    BALSAMIC_version : str = BALSAMIC.__version__
    config_creation_date : str = "1999-01-01 00:00"

    class Config:
        validate_all = True

    @validator("analysis_dir")
    def dirpath_always_abspath(cls, value):
        if not Path(value).exists():
            raise ValidationError(f'Path does not exist!')
        else:
            return str(Path(value).resolve())

    @validator("log")
    def parse_analysis_to_log_path(cls, value, values, **kwargs):
        return str(Path(values.get("analysis_dir")) / values.get("case_id") / "logs") + "/"

    @validator("fastq_path")
    def parse_analysis_to_fastq_path(cls, value, values, **kwargs):
        return str(Path(values.get("analysis_dir")) / values.get("case_id") / "analysis" / "fastq") + "/"

    @validator("script")
    def parse_analysis_to_script_path(cls, value, values, **kwargs):
        return str(Path(values.get("analysis_dir")) / values.get("case_id") / "scripts")  + "/"

    @validator("result")
    def parse_analysis_to_result_path(cls, value, values, **kwargs):
        return str(Path(values.get("analysis_dir")) / values.get("case_id") / "analysis")  + "/"

    @validator("benchmark")
    def parse_analysis_to_benchmark_path(cls, value, values, **kwargs):
        return str(Path(values.get("analysis_dir")) / values.get("case_id") / "benchmarks") + "/"

    @validator("dag")
    def parse_analysis_to_dag_path(cls, value, values, **kwargs):
        return str(Path(values.get("analysis_dir")) / values.get("case_id") / f'{values.get("case_id")}_BALSAMIC_{BALSAMIC.__version__}_graph.pdf')

    @validator("config_creation_date")
    def datetime_as_string(cls, value):
        return datetime.now().strftime("%Y-%m-%d %H:%M")


class SampleInstanceModel(BaseModel):
    """Holds attributes for samples used in analysis"""
    file_prefix : str
    sample_type : Literal["tumor", "normal"] = Field(alias="type")
    readpair_suffix : List[str] = ["1", "2"]
    


class BioinfoToolsModel(BaseModel):
    """Holds versions of current bioinformatic tools used in analysis"""
    tabix : str
    bcftools : str
    fastqc : str
    manta : str
    picard : str
    bwa : str
    strelka : str
    gatk : str
    samtools : str
    sambamba : str
    vardic : Optional[str]
    cutadapt : Optional[str]


class PanelModel(BaseModel):
    capture_kit : Optional[str]
    chrom : Optional[List[str]]

    @validator("capture_kit")
    def path_as_abspath_str(cls, value):
        try:
            return str(Path(value).resolve())
        except:
            return None


class BalsamicConfigModel(BaseModel):
    """Concatenates config """
    QC : QCModel
    vcf : VCFModel
    analysis : AnalysisModel
    samples : Dict[str, SampleInstanceModel]
    reference : Dict[str, Path]
    conda_env_yaml : FilePath = str(CONDA_ENV_YAML)
    rule_directory : DirectoryPath = str(RULE_DIRECTORY) + "/"
    singularity : FilePath
    bioinfo_tools : BioinfoToolsModel
    panel : Optional[PanelModel]

    @validator("reference")
    def abspath_as_str(cls, value):
        for k, v in value.items():
            value[k] = str(Path(v).resolve())
        return value

    @validator("singularity")
    def transform_path_to_dict(cls, value):
        return {"image" : str(Path(value).resolve())}



