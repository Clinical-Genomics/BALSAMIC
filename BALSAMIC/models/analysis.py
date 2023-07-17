import hashlib
import logging
import re
from datetime import datetime
from pathlib import Path
from typing import Optional, List, Dict, Any

from pydantic import BaseModel, validator, Field, AnyUrl
from pydantic.types import DirectoryPath, FilePath

from BALSAMIC import __version__ as balsamic_version

from BALSAMIC.constants.analysis import (
    Gender,
    AnalysisType,
    AnalysisWorkflow,
    SequencingType,
    SampleType,
    MutationOrigin,
    MutationType,
    WorkflowSolution,
)
from BALSAMIC.constants.cache import GenomeVersion
from BALSAMIC.constants.constants import FileType
from BALSAMIC.constants.metrics import VALID_OPS

LOG = logging.getLogger(__name__)


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
        pop_freq: VCFAttributes (optional); maximum gnomad allele frequency
        pop_freq_umi: VCFAttributes (optional); maximum gnomad_af for UMI workflow
        strand_reads: VCFAttributes (optional); minimum strand specific read counts
        qss: VCFAttributes (optional); minimum sum of base quality scores
        sor: VCFAttributes (optional); minimum symmetrical log-odds ratio
        swegen_snv_freq: VCFAttributes (optional); maximum swegen snv allele frequency
        swegen_sv_freq: VCFAttributes (optional); maximum swegen sv allele frequency
        loqusdb_clinical_snv_freq: VCFAttributes (optional); maximum loqusdb clinical snv allele frequency
        loqusdb_clinical_sv_freq: VCFAttributes (optional); maximum loqusdb clinical sv allele frequency
        varcaller_name: str (required); variant caller name
        filter_type: str (required); filter name for variant caller
        analysis_type: str (required); analysis type e.g. tumor_normal or tumor_only
        description: str (required); comment section for description
    """

    AD: Optional[VCFAttributes]
    AF_min: Optional[VCFAttributes]
    AF_max: Optional[VCFAttributes]
    MQ: Optional[VCFAttributes]
    DP: Optional[VCFAttributes]
    pop_freq: Optional[VCFAttributes]
    pop_freq_umi: Optional[VCFAttributes]
    strand_reads: Optional[VCFAttributes]
    qss: Optional[VCFAttributes]
    sor: Optional[VCFAttributes]
    swegen_snv_freq: Optional[VCFAttributes]
    swegen_sv_freq: Optional[VCFAttributes]
    loqusdb_clinical_snv_freq: Optional[VCFAttributes]
    loqusdb_clinical_sv_freq: Optional[VCFAttributes]
    varcaller_name: str
    filter_type: str
    analysis_type: str
    description: str


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
        n_base_limit : Field(str(int)); supports filtering by limiting the N base number

    Raises:
        ValueError:
            When the input in min_seq_length and umi_trim_length cannot
            be interpreted as integer and coerced to string

    """

    picard_rmdup: bool = False
    adapter: str = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
    quality_trim: bool = True
    adapter_trim: bool = False
    umi_trim: bool = False
    min_seq_length: int = 25
    umi_trim_length: int = 5
    n_base_limit: int = 50

    @validator("min_seq_length", "umi_trim_length", "n_base_limit")
    def coerce_int_as_str(cls, value):
        return str(value)

    class Config:
        validate_all = True


class VarcallerAttribute(BaseModel):
    """Holds variables for variant caller software
    Attributes:
        mutation: str of mutation class
        mutation_type: str of mutation type
        analysis_type: list of str for analysis types
        workflow_solution: list of str for workflows
        sequencing_type: list of str for workflows

    Raises:
        ValueError:
            When a variable other than [somatic, germline] is passed in mutation field
            When a variable other than [SNV, CNV, SV] is passed in mutation_type field

    """

    mutation: str
    mutation_type: str = Field(alias="type")
    analysis_type: Optional[list]
    sequencing_type: Optional[list]
    workflow_solution: Optional[list]

    @validator("workflow_solution", check_fields=False)
    def workflow_solution_literal(cls, value) -> str:
        """Validate workflow solution"""
        assert set(value).issubset(
            set(WorkflowSolution)
        ), f"{value} is not valid workflow solution."
        return value

    @validator("analysis_type", check_fields=False)
    def annotation_type_literal(cls, value) -> str:
        """Validate analysis types"""
        assert set(value).issubset(
            set(AnalysisType)
        ), f"{value} is not a valid analysis type."
        return value

    @validator("mutation", check_fields=False)
    def mutation_literal(cls, value) -> str:
        """Validate mutation class"""
        assert value in set(MutationOrigin), f"{value} is not a valid mutation type."
        return value

    @validator("mutation_type", check_fields=False)
    def mutation_type_literal(cls, value) -> str:
        """Validate mutation type"""
        assert value in set(MutationType), f"{value} is not not a valid mutation class"
        return value

    @validator("sequencing_type", check_fields=False)
    def sequencing_type_literal(cls, value) -> str:
        """Validate sequencing type"""
        assert set(value).issubset(
            set(SequencingType)
        ), f"{value} is not not a valid sequencing type."
        return value


class VCFModel(BaseModel):
    """Contains VCF config"""

    vardict: VarcallerAttribute
    tnscope: VarcallerAttribute
    dnascope: VarcallerAttribute
    tnscope_umi: VarcallerAttribute
    manta_germline: VarcallerAttribute
    manta: VarcallerAttribute
    dellysv: VarcallerAttribute
    cnvkit: VarcallerAttribute
    ascat: VarcallerAttribute
    dellycnv: VarcallerAttribute
    tiddit: VarcallerAttribute
    cnvpytor: VarcallerAttribute
    svdb: VarcallerAttribute


class AnalysisModel(BaseModel):
    """Pydantic model containing workflow variables

    Attributes:

        case_id : Field(required); string case identifier
        gender: Field(required); string case gender
        analysis_type : Field(required); string literal [single, paired, pon]
            single : if only tumor samples are provided
            paired : if both tumor and normal samples are provided
            pon : panel of normal analysis
        sequencing_type : Field(required); string literal [targeted, wgs]
            targeted : if capture kit was used to enrich specific genomic regions
            wgs : if whole genome sequencing was performed
        analysis_workflow: Field(required); string literal [balsamic, balsamic-qc, balsamic-umi]
            balsamic: execute balsamic workflow
            balsamic-qc: execute balsamic qc-only workflow
            balsamic-umi: execute balsamic along with UMIworkflow for panels
        analysis_dir : Field(required); existing path where to save files
        fastq_path : Field(optional); Path where fastq files will be stored
        script : Field(optional); Path where snakemake scripts will be stored
        log : Field(optional); Path where logs will be saved
        result : Field(optional); Path where BALSAMIC output will be stored
        benchmark : Field(optional); Path where benchmark report will be stored
        dag : Field(optional); Path where DAG graph of workflow will be stored
        BALSAMIC_version  : Field(optional); Current version of BALSAMIC
        config_creation_date  : Field(optional); Timestamp when config was created

    Raises:
        ValueError:
            When gender is set to any other than [female, male]
            When analysis_type is set to any value other than [single, paired, pon]
            When sequencing_type is set to any value other than [wgs, targeted]
            When analysis_workflow is set to any other than [balsamic, balsamic-qc, balsamic-umi]
    """

    case_id: str
    analysis_type: str
    gender: Optional[str]
    sequencing_type: str
    analysis_workflow: str
    analysis_dir: DirectoryPath
    fastq_path: Optional[DirectoryPath]
    script: Optional[DirectoryPath]
    log: Optional[DirectoryPath]
    result: Optional[DirectoryPath]
    benchmark: Optional[DirectoryPath]
    dag: Optional[FilePath]
    BALSAMIC_version: str = balsamic_version
    config_creation_date: Optional[str]

    class Config:
        validate_all = True

    @validator("analysis_type")
    def analysis_type_literal(cls, value) -> str:
        if value not in set(AnalysisType):
            raise ValueError(
                f"Provided analysis type ({value}) not supported in BALSAMIC!"
            )
        return value

    @validator("gender")
    def gender_literal(cls, value, values) -> Optional[str]:
        if value not in set(Gender) and values.get("analysis_type") != "pon":
            raise ValueError(
                f"Provided gender type ({value}) is not supported in BALSAMIC!"
            )
        return value

    @validator("sequencing_type")
    def sequencing_type_literal(cls, value) -> str:
        if value not in set(SequencingType):
            raise ValueError(
                f"Provided sequencing type ({value}) not supported in BALSAMIC!"
            )
        return value

    @validator("analysis_workflow", check_fields=True)
    def analysis_workflow_literal(cls, value) -> str:
        if value not in set(AnalysisWorkflow):
            raise ValueError(
                f"Provided analysis workflow ({value} not supported in BALSAMIC"
            )
        return value

    @validator("analysis_dir")
    def dirpath_always_abspath(cls, value) -> str:
        return Path(value).resolve().as_posix()

    @validator("fastq_path")
    def fastq_path_as_abspath(cls, value, values, **kwargs) -> str:
        return Path(value).resolve().as_posix()

    @validator("log")
    def parse_analysis_to_log_path(cls, value, values, **kwargs) -> str:
        return (
            Path(values.get("analysis_dir"), values.get("case_id"), "logs").as_posix()
            + "/"
        )

    @validator("script")
    def parse_analysis_to_script_path(cls, value, values, **kwargs) -> str:
        return (
            Path(
                values.get("analysis_dir"), values.get("case_id"), "scripts"
            ).as_posix()
            + "/"
        )

    @validator("result")
    def parse_analysis_to_result_path(cls, value, values, **kwargs) -> str:
        return Path(
            values.get("analysis_dir"), values.get("case_id"), "analysis"
        ).as_posix()

    @validator("benchmark")
    def parse_analysis_to_benchmark_path(cls, value, values, **kwargs) -> str:
        return (
            Path(
                values.get("analysis_dir"), values.get("case_id"), "benchmarks"
            ).as_posix()
            + "/"
        )

    @validator("dag")
    def parse_analysis_to_dag_path(cls, value, values, **kwargs) -> str:
        return (
            Path(
                values.get("analysis_dir"), values.get("case_id"), values.get("case_id")
            ).as_posix()
            + f"_BALSAMIC_{balsamic_version}_graph.pdf"
        )

    @validator("config_creation_date")
    def datetime_as_string(cls, value):
        return datetime.now().strftime("%Y-%m-%d %H:%M")


class AnalysisPonModel(AnalysisModel):
    """Pydantic model containing PON workflow variables

    Attributes:
        pon_version: Field(str); version of the PON generated file
    """

    pon_version: str

    @validator("pon_version")
    def validate_pon_version(cls, value):
        """Checks that the version matches the following syntax: v<int>"""

        match = re.fullmatch("^v[1-9]\d*$", value)
        if not match:
            raise ValueError(
                f"The provided version ({value}) does not follow the defined syntax (v<int>)"
            )

        return value


class SampleInstanceModel(BaseModel):
    """Holds attributes for samples used in analysis.

    Attributes:
        type: Field(str): sample type [tumor, normal]
    """

    type: str

    @validator("type")
    def sample_type_literal(cls, value):
        """Validate balsamic supported sample type."""
        if value not in set(SampleType):
            raise ValueError(
                f"The provided sample type ({value}) is not supported in BALSAMIC"
            )
        return value


class PanelModel(BaseModel):
    """Holds attributes of PANEL BED file if provided
    Attributes:
        capture_kit : Field(str(Path)); string representation of path to PANEL BED file
        chrom : Field(list(str)); list of chromosomes in PANEL BED
        pon_cnn: Field(optional); Path where PON reference .cnn file is stored

    Raises:
        ValueError:
            When capture_kit argument is set, but is not a valid path

    """

    capture_kit: Optional[FilePath]
    chrom: Optional[List[str]]
    pon_cnn: Optional[FilePath]

    @validator("capture_kit")
    def path_as_abspath_str(cls, value):
        return Path(value).resolve().as_posix()

    @validator("pon_cnn")
    def pon_abspath_as_str(cls, value):
        if value:
            return Path(value).resolve().as_posix()
        return None


class PonBalsamicConfigModel(BaseModel):
    """Summarizes config models in preparation for export

    Attributes:
        QC : Field(QCmodel); variables relevant for fastq preprocessing and QC
        reference : Field(Dict); dictionary containing paths to reference genome files
        panel : Field(PanelModel(optional)); variables relevant to PANEL BED if capture kit is used
        singularity : Field(Path); path to singularity container of BALSAMIC
        rule_directory : Field(Path(RULE_DIRECTORY)); path where snakemake rules can be found
        bioinfo_tools : Field(dict); dictionary of bioinformatics software and which conda/container they are in
        bioinfo_tools_version : Field(dict); dictionary of bioinformatics software and their versions used for the analysis
    """

    QC: QCModel
    analysis: AnalysisPonModel
    samples: Dict[str, SampleInstanceModel]
    reference: Dict[str, Path]
    singularity: DirectoryPath
    bioinfo_tools: dict
    bioinfo_tools_version: dict
    panel: Optional[PanelModel]

    @validator("reference")
    def abspath_as_str(cls, value):
        for k, v in value.items():
            value[k] = Path(v).resolve().as_posix()
        return value

    @validator("singularity")
    def transform_path_to_dict(cls, value):
        return {"image": Path(value).resolve().as_posix()}


class BalsamicConfigModel(BaseModel):
    """Summarizes config models in preparation for export

    Attributes:
        QC : Field(QCmodel); variables relevant for fastq preprocessing and QC
        vcf : Field(VCFmodel); variables relevant for variant calling pipeline
        samples : Field(Dict); dictionary containing samples submitted for analysis
        reference : Field(Dict); dictionary containing paths to reference genome files
        panel : Field(PanelModel(optional)); variables relevant to PANEL BED if capture kit is used
        bioinfo_tools : Field(dict); dictionary of bioinformatics software and which conda/container they are in
        bioinfo_tools_version : Field(dict); dictionary of bioinformatics software and their versions used for the analysis
        singularity : Field(Path); path to singularity container of BALSAMIC
        background_variants: Field(Path(optional)); path to BACKGROUND VARIANTS for UMI
        rule_directory : Field(Path(RULE_DIRECTORY)); path where snakemake rules can be found
    """

    QC: QCModel
    vcf: Optional[VCFModel]
    analysis: AnalysisModel
    samples: Dict[str, SampleInstanceModel]
    reference: Dict[str, Path]
    singularity: DirectoryPath
    background_variants: Optional[FilePath]
    bioinfo_tools: dict
    bioinfo_tools_version: dict
    panel: Optional[PanelModel]

    @validator("reference")
    def abspath_as_str(cls, value):
        for k, v in value.items():
            value[k] = Path(v).resolve().as_posix()
        return value

    @validator("singularity")
    def transform_path_to_dict(cls, value):
        return {"image": Path(value).resolve().as_posix()}

    @validator("background_variants")
    def fl_abspath_as_str(cls, value):
        if value:
            return Path(value).resolve().as_posix()
        return None


class UMIParamsCommon(BaseModel):
    """This class defines the common params settings used as constants across various rules in UMI workflow.

    Attributes:
        align_format: str (required); output alignment format. eg. 'BAM'
        align_header: str (required); header line appended to the aligned BAM output
        align_intbases: int; input bases in each batch regardless of threads, for reproducibility
        filter_tumor_af: float (required); settings to filter minimum allelic frequency
    """

    align_header: str
    align_intbases: int
    filter_tumor_af: float


class UMIParamsUMIextract(BaseModel):
    """This class defines the params settings used as constants in UMI workflow-rule umextract.

    Attributes:
        read_structure: str (required); settings to define UMI read structure
    """

    read_structure: str = "-d, 'rs1,rs2'"


class UMIParamsConsensuscall(BaseModel):
    """This class defines the params settings used as constants in UMI workflow-rule consensuscall.

    Attributes:
        align_format: str (required); output alignment format. eg. 'BAM'
            filter_minreads: str (required); settings to filter consensus tags based on group size
        tag: str; Logic UMI tag
    """

    align_format: str = "BAM"
    filter_minreads: str = "3,1,1"
    tag: str = "XR"


class UMIParamsTNscope(BaseModel):
    """This class defines the params settings used as constants in UMI workflow- rule tnscope.

    Attributes:
        algo: str; choice of sentieon varcall algorithm. eg. 'TNscope'
        disable_detect: str; disable variant detector. eg 'sv' or 'snv_indel'
        filter_tumor_af: float (required); minimum allelic frequency to detect
        min_tumorLOD: int (required); minimum tumor log odds in the final call of variants
        init_tumorLOD: float (required); minimum tumor log odds in the initial pass calling variants
        error_rate: int (required); allow error-rate to consider in calling
        prunefactor: int (required); pruning factor in the kmer graph
        padding: int(required); amount to pad bed interval regions
    """

    algo: str
    init_tumorLOD: float
    min_tumorLOD: int
    error_rate: int
    prunefactor: int
    padding: int
    disable_detect: str


class ParamsVardict(BaseModel):
    """This class defines the params settings used as constants in vardict rule.

    Attributes:
        allelic_frequency: float (required); minimum allelic frequency to detect
        max_pval: float (required); the maximum p-value. Vardict default: 0.05
        max_mm: float (required); the maximum mean mismatches allowed. Vardict default: 5.25
        column_info: str (required); set of vardict filters for passing final variants
    """

    allelic_frequency: float
    max_pval: float
    max_mm: float
    column_info: str


class ParamsCommon(BaseModel):
    """This class defines the common params settings used as constants across various rules in balsamic workflow.

    Attributes:
        pcr_model: str (required). PCR indel model used to weed out false positive indels. Eg: none- PCR free samples.
        align_header: str (required); header line appended to the aligned BAM output
        min_mapq: int (required); minimum mapping quality score. Eg: 20- probability of mapping random read at 99% accuracy
        picard_fixmate: str (required), fix read mate information in bam file
        picard_RG_normal: str (required); replace readgroups in normal bam file
        picard_RG_tumor: str (required); replace readgroups in tumor bam file
    """

    align_header: str
    pcr_model: str
    min_mapq: int
    picard_fixmate: str
    picard_RG_normal: str
    picard_RG_tumor: str


class ParamsVEP(BaseModel):
    """This class defines the params settings used as constants in vep rule.

    Attributes:
        vep_filters: str (required); set of choosen options for processing vep annotated vcf file
    """

    vep_filters: str


class BalsamicWorkflowConfig(BaseModel):
    """Defines set of rules in balsamic workflow

    Handles attributes for corresponding rules.

    Attributes:
        common: global params defined across all rules in balsamic workflow
        umicommon: global params defined across specific rules in UMI workflow
        vep: global params defined in the rule vep
        vardict: params defined in the rule vardict
        umiextract : params defined in the rule sentieon_umiextract
        umiconsensuscall: params defined in the rule sentieon_consensuscall
        tnscope_umi: params defined in the rule sentieon_tnscope_umi
    """

    common: ParamsCommon
    vardict: ParamsVardict
    vep: ParamsVEP
    umicommon: UMIParamsCommon
    umiextract: UMIParamsUMIextract
    umiconsensuscall: UMIParamsConsensuscall
    tnscope_umi: UMIParamsTNscope