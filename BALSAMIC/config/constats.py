import BALSAMIC
import sys
import pydantic

from pydantic import BaseModel, ValidationError, validator
from pydantic.types import DirectoryPath, FilePath
from typing import Optional, List, Dict
from typing_extensions import Literal
from pathlib import Path

CONDA_ENV_PATH = Path(sys.modules["BALSAMIC"].__file__).parent.resolve() / "conda"
CONDA_ENV_YAML = Path(sys.modules["BALSAMIC"].__file__).parent.resolve() / "config" / "balsamic_env.yaml"
RULE_DIRECTORY = Path(sys.modules["BALSAMIC"].__file__).parent.resolve()


class ConfigQCModel(BaseModel):
    picard_rmdup : bool = False
    adapter: str = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
    quality_trim: bool = True
    adapter_trim: bool = False
    umi_trim: bool = False
    min_seq_length: str = "25"
    umi_trim_length: str = "5"



class ConfigVCFModel(BaseModel):
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


class ConfigAnalysisModel(BaseModel):
    case_id: str 
    analysis_type: Literal["single", "paired"]
    analysis_dir: Path
    fastq_path: Path = "analysis_dir/case_id/analysis/fastq_path"
    script: Path = "analysis_dir/case_id/scripts"
    log: Path = "analysis_dir/case_id/logs"
    result: Path = "analysis_dir/case_id/analysis"
    benchmark: Path = "analysis_dir/case_id/benchmarks"
    BALSAMIC_version : str = BALSAMIC.__version__

    class Config:
        validate_all = True

    @validator("analysis_dir")
    def dir_always_abspath(cls, value):
        if not Path(value).exists():
            raise ValidationError(f'Path does not exist!')
        else:
            return str(Path(value).resolve())

    @validator("log")
    def parse_analysis_to_log_path(cls, value, values, **kwargs):
        return str(Path(values.get("analysis_dir")) / values.get("case_id") / "logs")

    @validator("fastq_path")
    def parse_analysis_to_fastq_path(cls, value, values, **kwargs):
        return str(Path(values.get("analysis_dir")) / values.get("case_id") / "analysis" / "fastq")

    @validator("script")
    def parse_analysis_to_script_path(cls, value, values, **kwargs):
        return str(Path(values.get("analysis_dir")) / values.get("case_id") / "scripts")

    @validator("result")
    def parse_analysis_to_result_path(cls, value, values, **kwargs):
        return str(Path(values.get("analysis_dir")) / values.get("case_id") / "analysis")

    @validator("benchmark")
    def parse_analysis_to_benchmark_path(cls, value, values, **kwargs):
        return str(Path(values.get("analysis_dir")) / values.get("case_id") / "benchmarks")



class ConfigSampleInstanceModel(BaseModel):
    file_prefix : str
    sample_type : Literal["tumor", "normal"]
    readpair_suffix : List[str] = ["1", "2"]
    
    class Config:
        fields = {"sample_type": "type"}




class ConfigSamplesModel(BaseModel):
    class Config:
        extra = True


class ConfigReferenceModel(BaseModel):
    pass



class ConfigBioinfoToolsModel(BaseModel):
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


class BalsamicConfigModel(BaseModel):
    conda_env_yaml : FilePath = CONDA_ENV_YAML
    rule_directory : DirectoryPath = RULE_DIRECTORY
    singularity : FilePath 
    QC : ConfigQCModel
    vcf : ConfigVCFModel
    analysis : ConfigAnalysisModel
    samples : ConfigSamplesModel
    reference : ConfigReferenceModel
    bioinfo_tools : ConfigBioinfoToolsModel




