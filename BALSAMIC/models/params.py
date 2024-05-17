"""Balsamic analysis parameters models."""
from typing import Optional

from pydantic import BaseModel, ConfigDict
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


class ParamsManta(BaseModel):
    """This class defines the params settings used as constants in Manta rule.

    Attributes:
        wgs_settings: str(required). parameters for Manta analysis for WGS
        tga_settings: str(required). parameters for Manta analysis for TGA
    """

    wgs_settings: str
    tga_settings: str

class ParamsVEP(BaseModel):
    """This class defines the params settings used as constants in vep rule.

    Attributes:
        vep_filters: str (required); set of choosen options for processing vep annotated vcf file
    """

    vep_filters: str


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

    model_config = ConfigDict(coerce_numbers_to_str=True)
    picard_rmdup: bool = False
    adapter: str = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
    quality_trim: bool = True
    adapter_trim: bool = False
    umi_trim: bool = False
    min_seq_length: str = "25"
    umi_trim_length: str = "5"
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


class UMIParamsTNscope(BaseModel):
    """This class defines the params settings used as constants in UMI workflow-rule tnscope.

    Attributes:
        algo: str; choice of sentieon varcall algorithm. eg. 'TNscope'
        disable_detect: str; disable variant detector. eg 'sv' or 'snv_indel'
        filter_tumor_af: float (required); minimum allelic frequency to detect
        min_tumorLOD: int (required); minimum tumor log odds in the final call of variants
        init_tumorLOD: float (required); minimum tumor log odds in the initial pass calling variants
        error_rate: int (required); allow error-rate to consider in calling
        prunefactor: int (required); pruning factor in the kmer graph
        padding: int(required); amount to pad bed interval regions
        pcr_model: str (required). PCR indel model used to weed out false positive indels. Eg: none- PCR free samples.
    """

    algo: str
    filter_tumor_af: float
    init_tumorLOD: float
    min_tumorLOD: int
    error_rate: int
    prunefactor: int
    padding: int
    disable_detect: str
    pcr_model: str

class TGAParamsTNscope(BaseModel):
    """This class defines the params settings used as constants in TGA workflow-rule tnscope.

    Attributes:
        algo: str; choice of sentieon varcall algorithm. eg. 'TNscope'
        filter_tumor_af: float (required); minimum allelic frequency to detect
        min_tumorLOD: int (required); minimum tumor log odds in the final call of variants
        init_tumorLOD: float (required); minimum tumor log odds in the initial pass calling variants
        error_rate: int (required); allow error-rate to consider in calling
        prunefactor: int (required); pruning factor in the kmer graph
        padding: int(required); amount to pad bed interval regions
        pcr_model: str (required). PCR indel model used to weed out false positive indels. Eg: none- PCR free samples.
    """

    algo: str
    filter_tumor_af: float
    init_tumorLOD: float
    min_tumorLOD: int
    error_rate: int
    prunefactor: int
    padding: int
    pcr_model: str
class BalsamicWorkflowConfig(BaseModel):
    """Defines set of rules in balsamic workflow

    Handles attributes for corresponding rules.

    Attributes:
        common: global params defined across all rules in balsamic workflow
        manta: params used in the manta rules
        umicommon: global params defined across specific rules in UMI workflow
        vep: global params defined in the rule vep
        umiextract : params defined in the rule sentieon_umiextract
        umiconsensuscall: params defined in the rule sentieon_consensuscall
        tnscope_umi: params defined in the rule sentieon_tnscope_umi

    Functions:
        - get_manta_settings: Return setting for manta rule
    """

    common: ParamsCommon
    manta: ParamsManta
    vep: ParamsVEP
    umicommon: UMIParamsCommon
    umiextract: UMIParamsUMIextract
    umiconsensuscall: UMIParamsConsensuscall
    tnscope_umi: UMIParamsTNscope
    tnscope_tga:

    def get_manta_settings(self, sequencing_type) -> str:
        """Return correct setting for manta rules depending on sequencing type."""
        if sequencing_type == SequencingType.WGS:
            return self.manta.wgs_settings
        return self.manta.tga_settings


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
        high_normal_tumor_af_frac: VCFAttributes (optional); maximum normal allele frequency / tumor allele frequency
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
        low_pr_sr_count: VCFAttributes (optional); minumum Manta variant read support
        varcaller_name: str (required); variant caller name
        filter_type: str (required); filter name for variant caller
        analysis_type: str (required); analysis type e.g. tumor_normal or tumor_only
        description: str (required); comment section for description
    """

    AD: Optional[VCFAttributes] = None
    AF_min: Optional[VCFAttributes] = None
    high_normal_tumor_af_frac: Optional[VCFAttributes] = None
    DP: Optional[VCFAttributes] = None
    pop_freq: Optional[VCFAttributes] = None
    pop_freq_umi: Optional[VCFAttributes] = None
    strand_reads: Optional[VCFAttributes] = None
    qss: Optional[VCFAttributes] = None
    sor: Optional[VCFAttributes] = None
    swegen_snv_freq: Optional[VCFAttributes] = None
    swegen_sv_freq: Optional[VCFAttributes] = None
    loqusdb_clinical_snv_freq: Optional[VCFAttributes] = None
    loqusdb_clinical_sv_freq: Optional[VCFAttributes] = None
    low_pr_sr_count: Optional[VCFAttributes] = None
    varcaller_name: str
    filter_type: str
    analysis_type: str
    description: str
