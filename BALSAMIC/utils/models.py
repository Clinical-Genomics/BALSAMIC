import hashlib
from datetime import datetime
from pathlib import Path
from typing import Optional, List, Dict

from pydantic import BaseModel, validator, Field, AnyUrl
from pydantic.types import DirectoryPath, FilePath

from BALSAMIC import __version__ as balsamic_version

from BALSAMIC.utils.constants import (
    BIOINFO_TOOL_ENV, ANALYSIS_TYPES, WORKFLOW_SOLUTION, MUTATION_CLASS,
    MUTATION_TYPE, VALID_GENOME_VER, VALID_REF_FORMAT)


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
        pop_freq: VCFAttributes (optional); maximum gnomad_af
        strand_reads: VCFAttributes (optional); minimum strand specific read counts
        qss: VCFAttributes (optional); minimum sum of base quality scores
        sor: VCFAttributes (optional); minimum symmetrical log-odds ratio
        varcaller_name: str (required); variant caller name
        filter_type: str (required); filter name for variant caller
        analysis_type: str (required); analysis type e.g. tumor_normal or tumor_only
        description: str (required); comment section for description
    """

    AD: VCFAttributes
    AF_min: Optional[VCFAttributes]
    AF_max: Optional[VCFAttributes]
    MQ: Optional[VCFAttributes]
    DP: Optional[VCFAttributes]
    pop_freq: VCFAttributes
    strand_reads: Optional[VCFAttributes]
    qss: Optional[VCFAttributes]
    sor: Optional[VCFAttributes]
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

    @validator("min_seq_length", "umi_trim_length")
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

    Raises:
        ValueError:
            When a variable other than [somatic, germline] is passed in mutation field
            When a variable other than [SNV, CNV, SV] is passed in mutation_type field
            
    """

    mutation: str
    mutation_type: str = Field(alias="type")
    analysis_type: Optional[list]
    workflow_solution: Optional[list]

    @validator("workflow_solution", check_fields=False)
    def workflow_solution_literal(cls, value) -> str:
        " Validate workflow solution "
        assert set(value).issubset(
            set(WORKFLOW_SOLUTION)), f"{value} is not valid workflow solution."
        return value

    @validator("analysis_type", check_fields=False)
    def annotation_type_literal(cls, value) -> str:
        " Validate analysis types "
        assert set(value).issubset(
            set(ANALYSIS_TYPES)), f"{value} is not a valid analysis type."
        return value

    @validator("mutation", check_fields=False)
    def mutation_literal(cls, value) -> str:
        " Validate mutation class "
        assert value in MUTATION_CLASS, f"{value} is not a valid mutation type."
        return value

    @validator("mutation_type", check_fields=False)
    def mutation_type_literal(cls, value) -> str:
        " Validate mutation type "
        assert value in MUTATION_TYPE, f"{value} is not not a valid mutation class"
        return value


class VCFModel(BaseModel):
    """Contains VCF config"""

    manta: VarcallerAttribute
    cnvkit: VarcallerAttribute
    vardict: VarcallerAttribute
    tnscope: VarcallerAttribute
    dnascope: VarcallerAttribute
    tnhaplotyper: VarcallerAttribute
    manta_germline: VarcallerAttribute
    haplotypecaller: VarcallerAttribute
    TNscope_umi: VarcallerAttribute
    delly: VarcallerAttribute

class AnalysisModel(BaseModel):
    """Pydantic model containing workflow variables

    Attributes:
              
        case_id : Field(required); string case identifier
        analysis_type : Field(required); string literal [single, paired, pon]
            single : if only tumor samples are provided
            paired : if both tumor and normal samples are provided
            pon : panel of normal analysis
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

    Raises:
        ValueError:
            When analysis_type is set to any value other than [single, paired, qc, pon]
            When sequencing_type is set to any value other than [wgs, targeted]
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
    BALSAMIC_version: str = balsamic_version
    config_creation_date: Optional[str]

    class Config:
        validate_all = True

    @validator("analysis_type")
    def analysis_type_literal(cls, value) -> str:
        balsamic_analysis_types = ANALYSIS_TYPES
        if value not in balsamic_analysis_types:
            raise ValueError(
                f"Provided analysis type ({value}) not supported in BALSAMIC!")
        return value

    @validator("sequencing_type")
    def sequencing_type_literal(cls, value) -> str:
        balsamic_sequencing_types = ["wgs", "targeted"]
        if value not in balsamic_sequencing_types:
            raise ValueError(
                f"Provided sequencing type ({value}) not supported in BALSAMIC!"
            )
        return value

    @validator("analysis_dir")
    def dirpath_always_abspath(cls, value) -> str:
        return Path(value).resolve().as_posix()

    @validator("log")
    def parse_analysis_to_log_path(cls, value, values, **kwargs) -> str:
        return Path(values.get("analysis_dir"), values.get("case_id"),
                    "logs").as_posix() + "/"

    @validator("fastq_path")
    def parse_analysis_to_fastq_path(cls, value, values, **kwargs) -> str:
        return (Path(values.get("analysis_dir"), values.get("case_id"),
                     "analysis", "fastq").as_posix() + "/")

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
        return (Path(values.get("analysis_dir"), values.get("case_id"),
                     "benchmarks").as_posix() + "/")

    @validator("dag")
    def parse_analysis_to_dag_path(cls, value, values, **kwargs) -> str:
        return Path(values.get("analysis_dir"), values.get("case_id"),
                    values.get("case_id")).as_posix(
                    ) + f'_BALSAMIC_{balsamic_version}_graph.pdf'

    @validator("config_creation_date")
    def datetime_as_string(cls, value):
        return datetime.now().strftime("%Y-%m-%d %H:%M")


class SampleInstanceModel(BaseModel):
    """Holds attributes for samples used in analysis
    
    Attributes:
        file_prefix : Field(str); basename of sample pair
        sample_type : Field(str; alias=type); type of sample [tumor, normal]
        sample_name : Field(str); Internal ID of sample to use in deliverables
        readpair_suffix : Field(List); currently always set to [1, 2]
    
    Raises:
        ValueError:
            When sample_type is set ot any value other than [tumor, normal]

        """

    file_prefix: str
    sample_name: Optional[str]
    sample_type: str = Field(alias="type")
    readpair_suffix: List[str] = ["1", "2"]

    @validator("sample_type")
    def sample_type_literal(cls, value):
        balsamic_sample_types = ["tumor", "normal"]
        if value not in balsamic_sample_types:
            raise ValueError(
                f"Provided sample type ({value}) not supported in BALSAMIC!")
        return value

    @validator("sample_name")
    def set_sample_id_if_missing_value(cls, value, values, **kwargs):
        if value:
            return value
        return values.get("file_prefix")


class PanelModel(BaseModel):
    """Holds attributes of PANEL BED file if provided
    Attributes:
        capture_kit : Field(str(Path)); string representation of path to PANEL BED file
        chrom : Field(list(str)); list of chromosomes in PANEL BED

    Raises:
        ValueError:
            When capture_kit argument is set, but is not a valid path

    """

    capture_kit: Optional[FilePath]
    chrom: Optional[List[str]]

    @validator("capture_kit")
    def path_as_abspath_str(cls, value):
        return Path(value).resolve().as_posix()


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
    analysis: AnalysisModel
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
        umiworkflow : Field(bool); whether UMI workflow to run in parallel
    """

    QC: QCModel
    vcf: VCFModel
    analysis: AnalysisModel
    samples: Dict[str, SampleInstanceModel]
    reference: Dict[str, Path]
    singularity: DirectoryPath
    background_variants: Optional[FilePath]
    bioinfo_tools: dict
    bioinfo_tools_version: dict
    panel: Optional[PanelModel]
    umiworkflow: bool

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


class ReferenceUrlsModel(BaseModel):
    """Defines a basemodel for reference urls
    
    This class handles four attributes for each reference url. Each attribute defines url, type of file, and gzip status.

    Attributes:
      url: defines the url to access file. Essentially it will be used to download file locally. It should match url_type://...
      file_type: describes file type. Accepted values are VALID_REF_FORMAT constant 
      gzip: gzip status. Binary: True or False
      genome_version: genome version matching the content of the file. Accepted values are VALID_GENOME_VER constant 

    Raises:
      ValidationError: When it can't validate values matching above attributes
      
    """

    url: AnyUrl
    file_type: str
    gzip: bool = True
    genome_version: str
    output_file: Optional[str]
    output_path: Optional[str]
    secret: Optional[str]

    @validator("file_type")
    def check_file_type(cls, value) -> str:
        """Validate file format according to constants"""
        assert value in VALID_REF_FORMAT, f"{value} not a valid reference file format."
        return value

    @validator("genome_version")
    def check_genome_ver(cls, value) -> str:
        """Validate genome version according constants"""
        assert value in VALID_GENOME_VER, f"{value} not a valid genome version."
        return value

    @property
    def get_output_file(self):
        """return output file full path"""
        output_file_path = Path(self.output_path, self.output_file).as_posix()
        return output_file_path

    @property
    def write_md5(self):
        """calculate md5 for first 4kb of file and write to file_name.md5"""
        hash_md5 = hashlib.md5()
        output_file = Path(self.output_path, self.output_file)
        if not output_file.is_file():
            raise FileNotFoundError(
                f"{output_file.as_posix()} file does not exist")

        with open(output_file.as_posix(), "rb") as fh:
            for chunk in iter(lambda: fh.read(4096), b""):
                hash_md5.update(chunk)

        with open(output_file.as_posix() + ".md5", "w") as fh:
            fh.write("{} {}\n".format(output_file.as_posix(),
                                      hash_md5.hexdigest()))


class ReferenceMeta(BaseModel):
    """Defines a basemodel for all reference file

    This class defines a meta for various reference files. Only reference_genome is mandatory.

    Attributes:
        basedir: str for base directory which will be appended to all ReferenceUrlsModel fields
        reference_genome: ReferenceUrlsModel. Required field for reference genome fasta file
        dbsnp: ReferenceUrlsModel. Optional field for dbSNP vcf file
        hc_vcf_1kg: ReferenceUrlsModel. Optional field for high confidence 1000Genome vcf
        mills_1kg: ReferenceUrlsModel. Optional field for Mills' high confidence indels vcf
        known_indel_1kg: ReferenceUrlsModel. Optional field for 1000Genome known indel vcf
        vcf_1kg: ReferenceUrlsModel. Optional field for 1000Genome all SNPs
        wgs_calling: ReferenceUrlsModel. Optional field for wgs calling intervals
        genome_chrom_size: ReferenceUrlsModel. Optional field for geneome's chromosome sizes
        gnomad_variant: ReferenceUrlsModel. Optional gnomad variants (non SV) as vcf
        cosmicdb: ReferenceUrlsModel. Optional COSMIC database's variants as vcf
        refgene_txt: ReferenceUrlsModel. Optional refseq's gene flat format from UCSC
        refgene_sql: ReferenceUrlsModel. Optional refseq's gene sql format from UCSC
        rankscore: ReferenceUrlsModel. Optional rankscore model
        access_regions: ReferenceUrlsModel. Optional field for accessible genome regions
        delly_exclusion: ReferenceUrlsModel. Optional field for genome exclusion regions
        ascat_gccorrection: ReferenceUrlsModel. Optional field for genome gc correction bins
        ascat_chryloci: ReferenceUrlsModel. Optional field for chromosome Y loci
    """

    basedir: str = ""
    reference_genome: ReferenceUrlsModel
    dbsnp: Optional[ReferenceUrlsModel]
    hc_vcf_1kg: Optional[ReferenceUrlsModel]
    mills_1kg: Optional[ReferenceUrlsModel]
    known_indel_1kg: Optional[ReferenceUrlsModel]
    vcf_1kg: Optional[ReferenceUrlsModel]
    wgs_calling: Optional[ReferenceUrlsModel]
    genome_chrom_size: Optional[ReferenceUrlsModel]
    gnomad_variant: Optional[ReferenceUrlsModel]
    gnomad_variant_index: Optional[ReferenceUrlsModel]
    cosmicdb: Optional[ReferenceUrlsModel]
    refgene_txt: Optional[ReferenceUrlsModel]
    refgene_sql: Optional[ReferenceUrlsModel]
    rankscore: Optional[ReferenceUrlsModel]
    access_regions: Optional[ReferenceUrlsModel]
    delly_exclusion: Optional[ReferenceUrlsModel]
    ascat_gccorrection: Optional[ReferenceUrlsModel]
    ascat_chryloci: Optional[ReferenceUrlsModel]

    @validator("*", pre=True)
    def validate_path(cls, value, values, **kwargs):
        """validate and append path in ReferenceUrlsModel fields with basedir"""
        if isinstance(value, str):
            output_value = value
        else:
            if "output_path" in value:
                value["output_path"] = Path(values["basedir"],
                                            value["output_path"]).as_posix()
                output_value = ReferenceUrlsModel.parse_obj(value)
            else:
                output_value = value

        return output_value


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

    align_format: str = 'BAM'
    filter_minreads: str = '3,1,1'
    tag: str = 'XR'


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
    """

    algo: str
    init_tumorLOD: float
    min_tumorLOD: int
    error_rate: int
    prunefactor: int
    disable_detect: str


class UMIParamsVardict(BaseModel):
    """This class defines the params settings used as constants in UMIworkflow-rule vardict.

    Attributes:
        vardict_filters: str (required); set of filters to apply for variant-calling using vardict
    """
    vardict_filters: str


class UMIParamsVEP(BaseModel):
    """This class defines the params settings used as constants in UMIworkflow-rule vep.

    Attributes:
        vep_filters: str (required); set of filters to apply for variant-calling using vardict
    """
    vep_filters: str


class UMIworkflowConfig(BaseModel):
    """ Defines set of rules in UMI workflow.
 
    Handles attributes for corresponding rules.
	
    Attributes:
	common: global params defined across all rules in UMI workflow
	umiextract: params defined in the rule sentieon_umiextract
	consensuscall: params defined in the rule sentieon_consensuscall
	tnscope: params defined in the rule sentieon_tnscope_umi
	vardict: params defined in the rule vardict_umi
	vep: params defined in the rule vep_umi
    """
    common: UMIParamsCommon
    umiextract: UMIParamsUMIextract
    consensuscall: UMIParamsConsensuscall
    tnscope: UMIParamsTNscope
    vardict: UMIParamsVardict
    vep: UMIParamsVEP
