"""Balsamic analysis config case models."""

import re
from glob import glob
from pathlib import Path
from typing import Annotated, Dict, List, Optional

from pydantic import AfterValidator, BaseModel, field_validator, model_validator

from BALSAMIC import __version__ as balsamic_version
from BALSAMIC.constants.analysis import (
    AnalysisType,
    AnalysisWorkflow,
    FastqName,
    Gender,
    MutationOrigin,
    MutationType,
    PONWorkflow,
    SampleType,
    SequencingType,
    WorkflowSolution,
)
from BALSAMIC.models.params import QCModel
from BALSAMIC.models.validators import is_dir, is_file


class FastqInfoModel(BaseModel):
    """Holds filepaths for forward and reverse reads for a fastq_pattern."""

    fwd: Annotated[str, AfterValidator(is_file)]
    rev: Annotated[str, AfterValidator(is_file)]
    fwd_resolved: Annotated[Optional[str], AfterValidator(is_file)] = None
    rev_resolved: Annotated[Optional[str], AfterValidator(is_file)] = None


class SampleInstanceModel(BaseModel):
    """Holds attributes for samples used in analysis.

    Attributes:
        type: Field(str): sample type [tumor, normal]
        name: Field(str): sample name
        fastq_info: Field(dict): fastq patterns: paths to forward and reverse fastqs
    """

    type: SampleType
    name: str
    fastq_info: Dict[str, FastqInfoModel]


class PanelModel(BaseModel):
    """Holds attributes of PANEL BED file if provided
    Attributes:
        exome: (bool); optional parameter for targeted analyses to use exome parameters
        capture_kit : Field(str(Path)); string representation of path to PANEL BED file
        chrom : Field(list(str)); list of chromosomes in PANEL BED
        pon_cnn: Field(optional); Path where PON reference .cnn file is stored

    Raises:
        ValueError:
            When capture_kit argument is set, but is not a valid path

    """

    exome: Optional[bool] = False
    capture_kit: Annotated[Optional[str], AfterValidator(is_file)] = None
    chrom: Optional[List[str]] = None
    pon_cnn: Annotated[Optional[str], AfterValidator(is_file)] = None


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

    mutation: MutationOrigin
    mutation_type: MutationType
    analysis_type: Optional[List[AnalysisType]] = None
    sequencing_type: Optional[List[SequencingType]] = None
    workflow_solution: Optional[List[WorkflowSolution]] = None


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
    igh_dux4: VarcallerAttribute
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
    analysis_type: AnalysisType
    gender: Optional[Gender] = None
    sequencing_type: SequencingType
    analysis_workflow: AnalysisWorkflow
    analysis_dir: Annotated[str, AfterValidator(is_dir)]
    fastq_path: Annotated[str, AfterValidator(is_dir)]
    log: Annotated[str, AfterValidator(is_dir)]
    script: Annotated[str, AfterValidator(is_dir)]
    result: Annotated[str, AfterValidator(is_dir)]
    benchmark: Annotated[str, AfterValidator(is_dir)]
    dag: str
    BALSAMIC_version: str = balsamic_version
    config_creation_date: str
    pon_version: Optional[str] = None
    pon_workflow: Optional[PONWorkflow] = None

    @field_validator("pon_version")
    def validate_pon_version(cls, pon_version: Optional[str]):
        """Checks that the PON version matches the following syntax: v<int>"""
        if pon_version and not re.fullmatch("^v[1-9]\d*$", pon_version):
            raise ValueError(
                f"The provided PON version ({pon_version}) does not follow the defined syntax (v<int>)"
            )
        return pon_version


class CustomFilters(BaseModel):
    """Variant calling custom filters."""

    umi_min_reads: str | None = None


class Sentieon(BaseModel):
    """
    Class providing common functions and variables for different balsamic workflows.

    Attributes:
        sentieon_install_dir: Field(required); path to Sentieon installation directory
        sentieon_exec:  Field(required); path to Sentieon executeable
        sentieon_license: Field(required); Sentieon license string
    """

    sentieon_install_dir: Annotated[str, AfterValidator(is_dir)]
    sentieon_exec: Annotated[str, AfterValidator(is_file)]
    sentieon_license: str
    dnascope_model: Annotated[str, AfterValidator(is_file)]
    tnscope_model: Annotated[str, AfterValidator(is_file)]


class ConfigModel(BaseModel):
    """
    Class providing common functions and variables for different balsamic workflows.

    Attributes:
        QC : Field(QCmodel); variables relevant for fastq preprocessing and QC
        samples : Field(List[SampleInstanceModel]); List containing samples submitted for analysis
        reference : Field(Dict); dictionary containing paths to reference genome files
        panel : Field(PanelModel(optional)); variables relevant to PANEL BED if capture kit is used
        bioinfo_tools : Field(dict); dictionary of bioinformatics software and which conda/container they are in
        bioinfo_tools_version : Field(dict); dictionary of bioinformatics software and their versions used for the analysis
        singularity : Field(Dict); path to singularity container of BALSAMIC
        vcf : Field(VCFmodel); variables relevant for variant calling pipeline
        background_variants: Field(Path(optional)); path to BACKGROUND VARIANTS for UMI
        analysis: Field(AnalysisModel); Pydantic model containing workflow variables
        custom_filters: Field(CustomFilters); custom parameters for variant filtering
        sentieon: Field(required); Sentieon model attributes

    This class also contains functions that help retrieve sample and file information,
    facilitating BALSAMIC run operations in Snakemake.

    Functions:
        - get_all_sample_names: Return all sample names in the analysis.
        - get_fastq_patterns_by_sample: Return all fastq patterns for given samples.
        - get_all_fastqs_for_sample: Return all fastqs for a sample.
        - get_fastq_by_fastq_pattern: Return fastq file path for requested fastq pattern and type.
        - get_sample_name_by_type: Return sample name for requested sample type.
        - get_sample_type_by_name: Return sample type for requested sample name.
        - get_bam_name_per_lane: Return list of bam file names for all fastq patterns of a sample.
        - get_final_bam_name: Return final bam name for downstream analysis.
    """

    QC: QCModel
    samples: List[SampleInstanceModel]
    reference: Dict[str, Path]
    singularity: Dict[str, str]
    bioinfo_tools: Dict
    bioinfo_tools_version: Dict
    panel: Optional[PanelModel] = None
    vcf: Optional[VCFModel] = None
    background_variants: Optional[str] = None
    analysis: AnalysisModel
    custom_filters: CustomFilters | None = None
    sentieon: Sentieon

    @field_validator("reference")
    def abspath_as_str(cls, reference: Dict[str, Path]):
        for k, v in reference.items():
            reference[k] = Path(v).resolve().as_posix()
        return reference

    @field_validator("singularity")
    def transform_path_to_dict(cls, singularity: Dict[str, str]):
        for k, v in singularity.items():
            singularity[k] = Path(v).resolve().as_posix()
        return singularity

    @field_validator("background_variants")
    def background_variants_abspath_as_str(cls, background_variants: str):
        """Converts FilePath to string."""
        if background_variants:
            return Path(background_variants).resolve().as_posix()
        return None

    @field_validator("samples")
    def no_duplicate_fastq_patterns(cls, samples):
        """Validate that no duplicate fastq patterns have been assigned in dict."""
        fastq_pattern_counts = {}

        # Count Fastq pattern occurrence
        for sample in samples:
            for fastq_pattern in sample.fastq_info.keys():
                if fastq_pattern not in fastq_pattern_counts:
                    fastq_pattern_counts[fastq_pattern] = 1
                else:
                    fastq_pattern_counts[fastq_pattern] += 1

        # Look for duplicates
        duplicates = []
        for fastq_pattern in fastq_pattern_counts:
            if fastq_pattern_counts[fastq_pattern] > 1:
                duplicates.append(fastq_pattern)

        if duplicates:
            raise ValueError(
                f"Duplicate FastqPattern(s) found: {', '.join(duplicates)} across multiple samples"
            )

        return samples

    @model_validator(mode="before")
    def no_unassigned_fastqs_in_fastq_dir(cls, values):
        """All fastq files in the supplied fastq-dir must have been assigned to the sample-dict."""

        def get_all_fwd_rev_values(samples) -> List[str]:
            # Return all fastq files in analysis
            fwd_rev_values = []
            for sample in samples:
                for fastq_pattern in sample["fastq_info"]:
                    fwd_rev_values.append(
                        sample["fastq_info"][fastq_pattern][FastqName.FWD]
                    )
                    fwd_rev_values.append(
                        sample["fastq_info"][fastq_pattern][FastqName.REV]
                    )
            return fwd_rev_values

        fastq_path = values["analysis"]["fastq_path"]

        # Get a set of all fastq files in fastq-directory
        fastqs_in_fastq_path = set(glob(f"{fastq_path}/*fastq.gz"))

        # Look for fastqs in sample dict
        fastqs_assigned = set(get_all_fwd_rev_values(values["samples"]))

        unassigned_fastqs = fastqs_in_fastq_path - fastqs_assigned
        if unassigned_fastqs:
            raise ValueError(
                f"Fastqs in fastq-dir not assigned to sample config: {unassigned_fastqs}"
            )

        return values

    def get_all_sample_names(self) -> List[str]:
        """Return all sample names in the analysis."""
        return [sample.name for sample in self.samples]

    def get_fastq_patterns_by_sample(self, sample_names: List[str]) -> List[str]:
        """Return all fastq_patterns for a given sample."""
        return [
            fastq_pattern
            for sample in self.samples
            if sample.name in sample_names
            for fastq_pattern in sample.fastq_info.keys()
        ]

    def get_all_fastqs_for_sample(
        self, sample_name: str, fastq_types: Optional[List[FastqName]] = None
    ) -> List[str]:
        """Return all fastqs (optionally only [fwd/rev]) involved in analysis of sample."""

        if fastq_types is None:
            fastq_types = [FastqName.FWD, FastqName.REV]

        fastq_list: List = []
        for sample in self.samples:
            if sample.name == sample_name:
                for fastq_info in sample.fastq_info.values():
                    if FastqName.FWD in fastq_types:
                        fastq_list.append(fastq_info.fwd)
                    if FastqName.REV in fastq_types:
                        fastq_list.append(fastq_info.rev)
        return fastq_list

    def get_all_fastq_names(self, remove_suffix: bool = False) -> List[str]:
        """Return all fastq_names involved in analysis, optionally remove fastq.gz suffix."""
        fastq_names = []
        for sample in self.samples:
            for fastq_pattern, fastqs in sample.fastq_info.items():
                if remove_suffix:
                    fastq_names.extend(
                        [
                            Path(fastqs.fwd).name.replace(".fastq.gz", ""),
                            Path(fastqs.rev).name.replace(".fastq.gz", ""),
                        ]
                    )
                else:
                    fastq_names.extend(
                        [
                            Path(fastqs.fwd).name,
                            Path(fastqs.rev).name,
                        ]
                    )
        return fastq_names

    def get_fastq_by_fastq_pattern(
        self, fastq_pattern: str, fastq_type: FastqName
    ) -> str:
        """Return fastq file path for requested fastq pair pattern and fastq type: [fwd/rev]."""
        for sample in self.samples:
            if fastq_pattern in sample.fastq_info:
                return (
                    sample.fastq_info[fastq_pattern].fwd
                    if fastq_type == FastqName.FWD
                    else sample.fastq_info[fastq_pattern].rev
                )

    def get_sample_name_by_type(self, sample_type: str) -> str:
        """Return sample name for requested sample type."""
        for sample in self.samples:
            if sample.type == sample_type:
                return sample.name

    def get_sample_type_by_name(self, sample_name: str, uppercase: bool = False) -> str:
        """Return sample type for requested sample name, optionally return it capitalized"""
        for sample in self.samples:
            if sample.name == sample_name:
                return sample.type.upper() if uppercase else sample.type

    def get_bam_name_per_lane(self, bam_dir: str, sample_name: str) -> List[str]:
        """Return list of bam-file names for all fastq_patterns of a sample."""
        bam_names = []
        for sample in self.samples:
            if sample.name == sample_name:
                bam_names.extend(
                    [
                        f"{bam_dir}{sample_name}_align_sort_{fastq_pattern}.bam"
                        for fastq_pattern in sample.fastq_info
                    ]
                )
        return bam_names

    def get_final_bam_name(
        self,
        bam_dir: str,
        sample_name: str = None,
        sample_type: str = None,
        specified_suffix: str = None,
    ) -> str:
        """Return final bam name to be used in downstream analysis."""

        if not sample_name and not sample_type:
            raise ValueError(
                "Either sample_name or sample_type must be provided to get the final bam name."
            )

        sample_name = (
            self.get_sample_name_by_type(sample_type)
            if not sample_name
            else sample_name
        )

        sample_type = (
            self.get_sample_type_by_name(sample_name)
            if not sample_type
            else sample_type
        )

        if self.analysis.analysis_type == AnalysisType.PON:
            # Only dedup is necessary for panel of normals
            final_bam_suffix = "dedup"
        elif self.analysis.sequencing_type == SequencingType.TARGETED:
            # Only dedup is necessary for TGA
            final_bam_suffix = "dedup.fixmate"
        else:
            # For WGS the bamfiles are realigned
            final_bam_suffix = "dedup.realign"

        if specified_suffix:
            final_bam_suffix = specified_suffix

        return f"{bam_dir}{sample_type}.{sample_name}.{final_bam_suffix}.bam"

    def get_cnv_report_plots(self) -> List[str]:
        """Return a list of AscatNgs CNV plot files."""
        if self.analysis.analysis_type == AnalysisType.SINGLE:
            return [
                f"CNV.somatic.{self.analysis.case_id}.cnvpytor.circular.png",
                f"CNV.somatic.{self.analysis.case_id}.cnvpytor.scatter.png",
            ]
        return [
            f"CNV.somatic.{self.analysis.case_id}.ascat.ascatprofile.png",
            f"CNV.somatic.{self.analysis.case_id}.ascat.rawprofile.png",
            f"CNV.somatic.{self.analysis.case_id}.ascat.ASPCF.png",
            f"CNV.somatic.{self.analysis.case_id}.ascat.tumor.png",
            f"CNV.somatic.{self.analysis.case_id}.ascat.germline.png",
            f"CNV.somatic.{self.analysis.case_id}.ascat.sunrise.png",
        ]
