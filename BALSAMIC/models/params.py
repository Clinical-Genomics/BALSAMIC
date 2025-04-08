"""Balsamic analysis parameters models."""
from typing import List, Optional

from pydantic import BaseModel, ConfigDict, field_validator
from BALSAMIC.constants.analysis import SequencingType


class ParamsCommon(BaseModel):
    """This class defines the common params settings used as constants across various rules in balsamic workflow.

    Attributes:
        pcr_model: str (required). PCR indel model used to weed out false positive indels. Eg: none- PCR free samples.
        min_mapq: int (required); minimum mapping quality score. Eg: 20- probability of mapping random read at 99% accuracy
        picard_fixmate: str (required), fix read mate information in bam file
        picard_RG_normal: str (required); replace readgroups in normal bam file
        picard_RG_tumor: str (required); replace readgroups in tumor bam file
    """

    pcr_model: str
    min_mapq: int
    picard_fixmate: str
    picard_RG_normal: str
    picard_RG_tumor: str


class ParamsInsertSizeMetrics(BaseModel):
    """This class defines the common params settings used for the InsertSizeMetricsAlgo

    Attributes:
        min_read_ratio: float(required). Minimum ratio of reads for a read category to be included in the output histogram
    """

    min_read_ratio: float


class ParamsManta(BaseModel):
    """This class defines the params settings used as constants in Manta rule.

    Attributes:
        wgs_settings: str(required). parameters for Manta analysis for WGS
        tga_settings: str(required). parameters for Manta analysis for TGA
    """

    wgs_settings: str
    tga_settings: str


class ParamsMosdepth(BaseModel):
    """This class defines the params settings used as constants in Mosdepth rule.

    Attributes:
        mapq: str(required); mapping quality threshold, reads with a quality less than this value are ignored
        samflag: str(required); exclude reads with any of the bits in FLAG set
        quantize: str(required); merges adjacent bases as long as they fall in the same coverage bins e.g. (10-20)
    """

    mapq: int
    samflag: int
    quantize: str


class ParamsSentieonWGSMetrics(BaseModel):
    """This class defines the params settings used as constants in Sentieon WGS Metrics rule.

    Attributes:
        min_base_qual: int(required); base quality threshold, bases with a quality less than this value are ignored
        cov_threshold: list(required); coverage threshold list
    """

    min_base_qual: int
    cov_threshold: str

    @field_validator("cov_threshold", mode="before")
    def parse_into_arguments(cls, cov_threshold):
        param_values = []
        for value in cov_threshold:
            param_values.append(f"--cov_thresh {value}")
        return " ".join(param_values)


class ParamsVEP(BaseModel):
    """This class defines the params settings used as constants in vep rule.

    Attributes:
        vep_filters: str (required); set of choosen options for processing vep annotated vcf file
    """

    vep_filters: str


class QCModel(BaseModel):
    """Contains settings for quality control and pre-processing
    Attributes:
        adapter : Field(str(AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT)); adapter sequence to trim
        min_seq_length : Field(str(int)); minimum sequence length cutoff for reads
        n_base_limit : Field(str(int)); supports filtering by limiting the N base number

    Raises:
        ValueError:
            When the input in min_seq_length and umi_trim_length cannot
            be interpreted as integer and coerced to string

    """

    model_config = ConfigDict(coerce_numbers_to_str=True)
    adapter: str = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
    min_seq_length: str = "25"
    n_base_limit: str = "50"


class UMIParamsCommon(BaseModel):
    """This class defines the common params settings used as constants across various rules in UMI workflow.

    Attributes:
        align_intbases: int; input bases in each batch regardless of threads, for reproducibility
    """

    align_intbases: int


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


class ParamsTNscope(BaseModel):
    """This class defines the TNscope params."""

    algo: str
    filter_tumor_af: Optional[float] = None
    init_tumorLOD: Optional[float] = None
    min_tumorLOD: Optional[float] = None
    init_normalLOD: Optional[float] = None
    min_normalLOD: Optional[float] = None
    error_rate: Optional[int] = None
    prunefactor: Optional[int] = None
    padding: Optional[int] = None
    disable_detect: Optional[str] = None
    pcr_model: Optional[str] = None


class BAMPostProcessingParams(BaseModel):
    """This class defines the params settings used as constants bam post-processing rules.

    Attributes:
       manta_max_base_quality: int (required); the maximum base quality in bamfile used downstream in Manta rules
    """

    manta_max_base_quality: int


class BEDPreProcessingParams(BaseModel):
    """This class defines the params settings used as constants in bed pre-processing rules.

    Attributes:
       minimum_region_size: int (required); the minimum region size in input bedfiles for CNV analysis
    """

    minimum_region_size: int


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


class BalsamicWorkflowConfig(BaseModel):
    """Defines set of rules in balsamic workflow

    Handles attributes for corresponding rules.

    Attributes:
        common: global params defined across all rules in balsamic workflow
        bam_post_processing: params used in bam post-processing rules
        manta: params used in the manta rules
        mosdepth: params used in mosdepth rule
        umicommon: global params defined across specific rules in UMI workflow
        vep: global params defined in the rule vep
        umiextract : params defined in the rule sentieon_umiextract
        umiconsensuscall: params defined in the rule sentieon_consensuscall
        tnscope_umi: params defined in the rule sentieon_tnscope_umi
        tnscope_tga_tumor_only: params defined in the rule sentieon_tnscope_tga_tumor_only
        tnscope_tga_tumor_normal: params defined in the rule sentieon_tnscope_tga_tumor_normal

    Functions:
        - get_manta_settings: Return setting for manta rule
    """

    bam_post_processing: BAMPostProcessingParams
    bed_pre_processing: BEDPreProcessingParams
    common: ParamsCommon
    insert_size_metrics: ParamsInsertSizeMetrics
    manta: ParamsManta
    mosdepth: ParamsMosdepth
    sentieon_wgs_metrics: ParamsSentieonWGSMetrics
    vep: ParamsVEP
    vardict: ParamsVardict
    umicommon: UMIParamsCommon
    umiextract: UMIParamsUMIextract
    umiconsensuscall: UMIParamsConsensuscall
    tnscope_umi: ParamsTNscope
    tnscope_tga: ParamsTNscope
    tnscope_wgs: ParamsTNscope

    def get_manta_settings(self, sequencing_type) -> str:
        """Return correct setting for manta rules depending on sequencing type."""
        if sequencing_type == SequencingType.WGS:
            return self.manta.wgs_settings
        return self.manta.tga_settings


class VCFFilter(BaseModel):
    """General purpose filter to manage various VCF attributes

    This class handles three parameters for the purpose filtering variants
    based on a tag_values, filter_name, and which field in the VCF.

    E.g. AD=VCFAttributes(tag_value=5, filter_name="balsamic_low_tumor_ad", field="INFO")
    A value of 5 from INFO field and filter_name will be balsamic_low_tumor_ad

    Attributes:
        tag_value: float
        filter_name: str
        field: str
        Description: str (optional); filter description
        analysis_type: str (optional); specific sequencing type such paired or single for which to apply the filter
        variant_caller: str (optional); the specific variant caller for which to apply the filter
        exome: bool (optional); if the filter should only be applied for exome
    """

    tag_value: Optional[float] = None
    filter_name: str
    field: Optional[str] = None
    Description: Optional[str] = None
    analysis_type: Optional[str] = None
    variant_caller: Optional[str] = None
    exome: Optional[bool] = None


class StructuralVariantFilters(BaseModel):
    """Variant filters for Structural Variants

    This class handles attributes and filter for structural variants

    Attributes:
        swegen_sv_freq: VCFAttributes (optional); maximum swegen sv allele frequency
        loqusdb_clinical_sv_freq: VCFAttributes (optional); maximum loqusdb clinical sv allele frequency
        low_pr_sr_count: VCFAttributes (optional); minumum Manta variant read support
        varcaller_name: str (required); variant caller name
        filter_type: str (required); filter name for variant caller
        analysis_type: str (required); analysis type e.g. tumor_normal or tumor_only
        description: str (required); comment section for description
    """

    swegen_sv_freq: Optional[VCFFilter] = None
    loqusdb_clinical_sv_freq: Optional[VCFFilter] = None
    low_pr_sr_count: Optional[VCFFilter] = None
    varcaller_name: str
    filter_type: str
    analysis_type: str
    description: str
