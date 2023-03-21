import hashlib
import logging
import re
import os
from datetime import datetime
from pathlib import Path
from typing import Optional, List, Dict, Any

from pydantic import BaseModel, validator, Field, AnyUrl, root_validator
from pydantic.types import DirectoryPath, FilePath

from BALSAMIC import __version__ as balsamic_version

from BALSAMIC.constants.common import (
    SEQUENCING_TYPE,
    ANALYSIS_TYPES,
    ANALYSIS_WORKFLOW,
    WORKFLOW_SOLUTION,
    MUTATION_CLASS,
    MUTATION_TYPE,
    VALID_OPS,
    GENDER_OPTIONS,
    SAMPLE_TYPE,
    FASTQ_SUFFIXES,
)
from BALSAMIC.constants.reference import VALID_GENOME_VER, VALID_REF_FORMAT

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
            set(WORKFLOW_SOLUTION)
        ), f"{value} is not valid workflow solution."
        return value

    @validator("analysis_type", check_fields=False)
    def annotation_type_literal(cls, value) -> str:
        """Validate analysis types"""
        assert set(value).issubset(
            set(ANALYSIS_TYPES)
        ), f"{value} is not a valid analysis type."
        return value

    @validator("mutation", check_fields=False)
    def mutation_literal(cls, value) -> str:
        """Validate mutation class"""
        assert value in MUTATION_CLASS, f"{value} is not a valid mutation type."
        return value

    @validator("mutation_type", check_fields=False)
    def mutation_type_literal(cls, value) -> str:
        """Validate mutation type"""
        assert value in MUTATION_TYPE, f"{value} is not not a valid mutation class"
        return value

    @validator("sequencing_type", check_fields=False)
    def sequencing_type_literal(cls, value) -> str:
        """Validate sequencing type"""
        assert set(value).issubset(
            set(SEQUENCING_TYPE)
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
        balsamic_analysis_types = ANALYSIS_TYPES
        if value not in balsamic_analysis_types:
            raise ValueError(
                f"Provided analysis type ({value}) not supported in BALSAMIC!"
            )
        return value

    @validator("gender")
    def gender_literal(cls, value, values) -> Optional[str]:
        if value not in GENDER_OPTIONS and values.get("analysis_type") != "pon":
            raise ValueError(
                f"Provided gender type ({value}) is not supported in BALSAMIC!"
            )
        return value

    @validator("sequencing_type")
    def sequencing_type_literal(cls, value) -> str:
        balsamic_sequencing_types = SEQUENCING_TYPE
        if value not in balsamic_sequencing_types:
            raise ValueError(
                f"Provided sequencing type ({value}) not supported in BALSAMIC!"
            )
        return value

    @validator("analysis_workflow", check_fields=True)
    def analysis_workflow_literal(cls, value) -> str:
        balsamic_analysis_workflow = ANALYSIS_WORKFLOW
        if value not in balsamic_analysis_workflow:
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
class FastqInfoModel(BaseModel):
    """Verifies that fastq files have been correctly assigned in fastq-pattern of sample.

    Requirements to pass fastq-inputs for a sample:
    - All fastq-files must exist
    - All forward and reverse fastq pairs must be of the same length
    - All forward and reverse fastq pairs must have exactly 1 difference in their names

    """
    fwd: FilePath
    rev: FilePath
    @root_validator(pre=False)
    def validate_fastq_pairs_in_fastq_pattern(cls, values):
        values["fwd"] = Path(values["fwd"]).resolve().as_posix()
        values["rev"] = Path(values["rev"]).resolve().as_posix()

        fwd_read_path = values["fwd"]
        rev_read_path = values["rev"]

        # Fastq files in a pair must have file-names of same length
        fwd_fastq_name = os.path.basename(fwd_read_path)
        rev_fastq_name = os.path.basename(rev_read_path)

        len_fwd_fastq_name = len(fwd_fastq_name)
        len_rev_fastq_name = len(rev_fastq_name)

        assert len_fwd_fastq_name == len_rev_fastq_name, f"Fastq pair does not have names of equal length: " \
                                                         f"{fwd_fastq_name} {rev_fastq_name}"

        # Fastq files in a pair must have file-names with maximum 1 character difference
        count_str_diff = 0
        for i in range(len_fwd_fastq_name):
            if fwd_fastq_name[i] != rev_fastq_name[i]:
                count_str_diff += 1

        assert count_str_diff == 1, f"Fastq pair does not have exactly 1 differences ({count_str_diff})," \
                                    f"Fwd: {fwd_fastq_name} Rev: {rev_fastq_name}"

        return values

class SampleInstanceModel(BaseModel):
    """Holds attributes for samples used in analysis.

    Attributes:
        type: Field(str): sample type [tumor, normal]
        fastq_info: Field(dict): fastq patterns: paths to forward and reverse fastqs
    """

    type: str
    fastq_info: Dict[str, FastqInfoModel]

    @validator("type")
    def sample_type_literal(cls, value):
        """Validate balsamic supported sample type."""
        if value not in SAMPLE_TYPE:
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
        samples : Field(Dict); dictionary containing samples and their respective fastq-inputs submitted for analysis
        reference : Field(Dict); dictionary containing paths to reference genome files
        panel : Field(PanelModel(optional)); variables relevant to PANEL BED if capture kit is used
        bioinfo_tools : Field(dict); dictionary of bioinformatics software and which conda/container they are in
        bioinfo_tools_version : Field(dict); dictionary of bioinformatics software and their versions used for the analysis
        singularity : Field(Path); path to singularity container of BALSAMIC
        background_variants: Field(Path(optional)); path to BACKGROUND VARIANTS for UMI
        rule_directory : Field(Path(RULE_DIRECTORY)); path where snakemake rules can be found

    Requirements to pass fastq-inputs for a group of samples:
        - Fastq pattern can only be assigned once per group of samples
        - Fastq name can only be assigned once per group of samples
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

    @validator("samples", pre=True)
    def verify_unique_fastq_info_across_samples(cls, value):
        """Cool docstring."""
        fastq_patterns = {}
        fastq_filenames = {}
        for sample in value:
            for fastq_pattern in value[sample]["fastq_info"]:
                # Count occurrences of fastq pattern
                if fastq_pattern in fastq_patterns:
                    fastq_patterns[fastq_pattern] += 1
                else:
                    fastq_patterns[fastq_pattern] = 1

                # Count occurrences of fastq names
                for fastq_read_direction, fastq_path in value[sample]["fastq_info"][fastq_pattern].items():
                    fastq_name = os.path.basename(fastq_path)
                    if fastq_name in fastq_filenames:
                        fastq_filenames[fastq_name] += 1
                    else:
                        fastq_filenames[fastq_name] = 1

        # Fastq pattern can only be assigned once per group of samples
        for fastq_pattern in fastq_patterns:
            assert fastq_patterns[fastq_pattern] == 1, f"Fastq-pattern: {fastq_pattern} has been assigned more than once."

        # Fastq name can only be assigned once per group of samples
        for fastq_name in fastq_filenames:
            assert fastq_filenames[fastq_name] == 1, f"Fastq-name: {fastq_name} has been assigned more than once."

        return value


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
            raise FileNotFoundError(f"{output_file.as_posix()} file does not exist")

        with open(output_file.as_posix(), "rb") as fh:
            for chunk in iter(lambda: fh.read(4096), b""):
                hash_md5.update(chunk)

        with open(output_file.as_posix() + ".md5", "w") as fh:
            fh.write("{} {}\n".format(output_file.as_posix(), hash_md5.hexdigest()))


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
        delly_mappability: ReferenceUrlsModel. Optional field for genome mappability
        ascat_gccorrection: ReferenceUrlsModel. Optional field for genome gc correction bins
        ascat_chryloci: ReferenceUrlsModel. Optional field for chromosome Y loci
        clinvar: ReferenceUrlsModel. Optional field for clinvar reference
        somalier_sites: ReferenceUrlsModel. Optional field for somalier sites vcf
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
    delly_mappability: Optional[ReferenceUrlsModel]
    delly_mappability_gindex: Optional[ReferenceUrlsModel]
    delly_mappability_findex: Optional[ReferenceUrlsModel]
    ascat_gccorrection: Optional[ReferenceUrlsModel]
    ascat_chryloci: Optional[ReferenceUrlsModel]
    clinvar: Optional[ReferenceUrlsModel]
    somalier_sites: Optional[ReferenceUrlsModel]

    @validator("*", pre=True)
    def validate_path(cls, value, values, **kwargs):
        """validate and append path in ReferenceUrlsModel fields with basedir"""
        if isinstance(value, str):
            output_value = value
        else:
            if "output_path" in value:
                value["output_path"] = Path(
                    values["basedir"], value["output_path"]
                ).as_posix()
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


class MetricConditionModel(BaseModel):
    """Defines the metric condition model

    Attributes:
        norm: string (optional); validation condition
        threshold: float (optional); validation cut off
    """

    norm: Optional[str] = None
    threshold: Optional[float] = None


class MetricModel(BaseModel):
    """Defines the metric attributes model

    Attributes:
        header: str (optional); data
        id: str (required); unique sample identifier (sample_id, case_id or project_id)
        input: str (required); input file
        name: str (required); metric name
        step: str (required); step that generated the metric
        value: Any (required and can take None as a value); metric value
        condition: MetricConditionModel (required and can take None as a value); metric validation condition
    """

    header: Optional[str]
    id: str
    input: str
    name: str
    step: str
    value: Any = ...
    condition: Optional[MetricConditionModel] = ...

    @validator("name")
    def validate_name(cls, name, values):
        """Updates the name if the source is FastQC"""

        if "fastqc-percent_duplicates" in name:
            return "PERCENT_DUPLICATION_R" + values["input"].split("_")[-2]

        return name


class MetricValidationModel(BaseModel):
    """Defines the metric validation model

    Attributes:
        metrics: List[MetricModel] (required); metric model to validate

    Raises:
        ValueError: when a metric does not meet its validation requirements
    """

    metrics: List[MetricModel]

    @validator("metrics", each_item=True)
    def validate_metrics(cls, metric):
        """Checks if a metric meets its filtering condition"""

        if metric.condition and not VALID_OPS[metric.condition.norm](
            metric.value, metric.condition.threshold
        ):
            raise ValueError(
                f"QC metric {metric.name}: {metric.value} validation has failed. "
                f"(Condition: {metric.condition.norm} {metric.condition.threshold}, ID: {metric.id})."
            )

        LOG.info(f"QC metric {metric.name}: {metric.value} meets its condition.")

        return metric
