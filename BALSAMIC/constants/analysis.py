"""Balsamic analysis workflow constants."""
from enum import StrEnum
from typing import Dict, List

from BALSAMIC.constants.cache import DockerContainers


class LogFile:
    """Logfile constants"""

    LOGNAME: str = "balsamic.log"


class SubmitSnakemake:
    """Constants for sbatch script running snakemake on cluster"""

    MAX_RUN_HOURS: int = 168


class SnakemakeDAG:
    """Constants for Snakemake DAG parsing and rendering."""

    DIGRAPH_HEADER: str = "digraph snakemake_dag {"
    HEADER: str = "snakemake_dag {"

    GRAPH_NAME: str = "BALSAMIC"
    GRAPH_LABEL_LOC: str = "t"

    GRAPHVIZ_FORMAT: str = "pdf"
    GRAPHVIZ_ENGINE: str = "dot"


class RunMode(StrEnum):
    """Balsamic workflow run mode."""

    CLUSTER: str = "cluster"
    LOCAL: str = "local"


RUN_MODES: List[RunMode] = [mode for mode in RunMode]


class Gender(StrEnum):
    """Gender options."""

    FEMALE: str = "female"
    MALE: str = "male"
    UNKNOWN: str = "unknown"


class AnalysisType(StrEnum):
    """Supported analysis types."""

    PAIRED: str = "paired"
    PON: str = "pon"
    SINGLE: str = "single"


class AnalysisWorkflow(StrEnum):
    """Available Balsamic workflows."""

    BALSAMIC: str = "balsamic"
    BALSAMIC_QC: str = "balsamic-qc"
    BALSAMIC_UMI: str = "balsamic-umi"


ANALYSIS_WORKFLOWS: List[AnalysisWorkflow] = [workflow for workflow in AnalysisWorkflow]


class AnnotationCategory(StrEnum):
    CLINICAL: str = "clinical"
    RESEARCH: str = "research"


VARIANT_OBSERVATION_METAVALUES = {
    "gnomad_variant": {
        "fields": ["AF", "AF_popmax"],
        "ops": ["self", "self"],
        "names": ["GNOMADAF", "GNOMADAF_popmax"],
        "category": AnnotationCategory.RESEARCH,
    },
    "clinvar": {
        "fields": [
            "CLNVID",
            "CLNREVSTAT",
            "CLNSIG",
            "ORIGIN",
            "ONC",
            "ONCDN",
            "ONCREVSTAT",
            "ONCDISDB",
            "ONCCONF",
            "CLNVC",
            "CLNVCSO",
        ],
        "ops": [
            "self",
            "self",
            "self",
            "self",
            "self",
            "self",
            "self",
            "self",
            "self",
            "self",
            "self",
        ],
        "names": [
            "CLNVID",
            "CLNREVSTAT",
            "CLNSIG",
            "ORIGIN",
            "ONC",
            "ONCDN",
            "ONCREVSTAT",
            "ONCDISDB",
            "ONCCONF",
            "CLNVC",
            "CLNVCSO",
        ],
        "category": AnnotationCategory.RESEARCH,
    },
    "cadd_snv": {
        "names": ["CADD"],
        "ops": ["mean"],
        "columns": [6],
        "category": AnnotationCategory.RESEARCH,
    },
    "swegen_snv_frequency": {
        "fields": ["AF", "AC_Hom", "AC_Het", "AC_Hemi"],
        "ops": ["self", "self", "self", "self"],
        "names": [
            "SWEGENAF",
            "SWEGENAAC_Hom",
            "SWEGENAAC_Het",
            "SWEGENAAC_Hemi",
        ],
        "category": AnnotationCategory.RESEARCH,
    },
    "artefact_snv_observations": {
        "fields": ["Frq", "Obs", "Hom"],
        "ops": ["self", "self", "self"],
        "names": ["ArtefactFrq", "ArtefactObs", "ArtefactHom"],
        "category": AnnotationCategory.CLINICAL,
    },
    "artefact_sv_observations": {
        "fields": ["Frq", "Obs", "Hom"],
        "ops": ["self", "self", "self"],
        "names": ["ArtefactFrq", "ArtefactObs", "ArtefactHom"],
        "category": AnnotationCategory.CLINICAL,
    },
    "clinical_snv_observations": {
        "fields": ["Frq", "Obs", "Hom"],
        "ops": ["self", "self", "self"],
        "names": ["Frq", "Obs", "Hom"],
        "category": AnnotationCategory.CLINICAL,
    },
    "cancer_germline_snv_observations": {
        "fields": ["Frq", "Obs", "Hom"],
        "ops": ["self", "self", "self"],
        "names": [
            "Cancer_Germline_Frq",
            "Cancer_Germline_Obs",
            "Cancer_Germline_Hom",
        ],
        "category": AnnotationCategory.CLINICAL,
    },
    "cancer_somatic_snv_observations": {
        "fields": ["Frq", "Obs", "Hom"],
        "ops": ["self", "self", "self"],
        "names": [
            "Cancer_Somatic_Frq",
            "Cancer_Somatic_Obs",
            "Cancer_Somatic_Hom",
        ],
        "category": AnnotationCategory.CLINICAL,
    },
    "cancer_somatic_snv_panel_observations": {
        "fields": ["Frq", "Obs", "Hom"],
        "ops": ["self", "self", "self"],
        "names": [
            "Cancer_Somatic_Panel_Frq",
            "Cancer_Somatic_Panel_Obs",
            "Cancer_Somatic_Panel_Hom",
        ],
        "category": AnnotationCategory.CLINICAL,
    },
}


class SequencingType(StrEnum):
    """Sequencing carried out."""

    TARGETED: str = "targeted"
    WGS: str = "wgs"
    WES: str = "wes"


class SampleType(StrEnum):
    """Balsamic sample type inputs."""

    NORMAL: str = "normal"
    TUMOR: str = "tumor"


class MutationOrigin(StrEnum):
    """Variations present in a sample."""

    GERMLINE: str = "germline"
    SOMATIC: str = "somatic"


class MutationType(StrEnum):
    """Types of variations present in a sample."""

    CNV: str = "CNV"
    SNV: str = "SNV"
    SV: str = "SV"


class WorkflowSolution(StrEnum):
    """Solution applied to a specific part of the analysis."""

    BALSAMIC: str = "BALSAMIC"
    DRAGEN: str = "DRAGEN"
    SENTIEON: str = "Sentieon"
    SENTIEON_UMI: str = "Sentieon_umi"


class PONWorkflow(StrEnum):
    """Panel Of Normal creation workflow type."""

    CNVKIT: str = "CNVkit"
    GENS_MALE: str = "GENS_male"
    GENS_FEMALE: str = "GENS_female"


PON_WORKFLOWS: List[PONWorkflow] = [workflow for workflow in PONWorkflow]


class BioinfoTools(StrEnum):
    """List of bioinformatics tools in Balsamic."""

    ASCAT: str = "ascatNgs"
    BCFTOOLS: str = "bcftools"
    BEDTOOLS: str = "bedtools"
    BGZIP: str = "bgzip"
    BWA: str = "bwa"
    CNVKIT: str = "cnvkit"
    CNVPYTOR: str = "cnvpytor"
    COMPRESS: str = "compress"
    CSVKIT: str = "csvkit"
    DELLY: str = "delly"
    VEP: str = "ensembl-vep"
    FASTP: str = "fastp"
    FASTQC: str = "fastqc"
    GATK: str = "gatk"
    GENMOD: str = "genmod"
    MANTA: str = "manta"
    MSISENSORPRO: str = "msisensorpro"
    MOSDEPTH: str = "mosdepth"
    MULTIQC: str = "multiqc"
    PICARD: str = "picard"
    SAMBAMBA: str = "sambamba"
    SAMTOOLS: str = "samtools"
    SOMALIER: str = "somalier"
    SVDB: str = "svdb"
    TABIX: str = "tabix"
    TIDDIT: str = "tiddit"
    TNSCOPE: str = "tnscope"
    VARDICT: str = "vardict"
    VCF2CYTOSURE: str = "vcf2cytosure"
    VCFANNO: str = "vcfanno"
    CADD: str = "cadd"
    PURECN: str = "purecn"
    PYSAM: str = "pysam"


class FastqName(StrEnum):
    """Fastq name parameters."""

    FWD: str = "fwd"
    REV: str = "rev"


FASTQ_SUFFIXES: Dict[str, Dict] = {
    "1": {"fwd": "_1.fastq.gz", "rev": "_2.fastq.gz"},
    "2": {"fwd": "_R1_001.fastq.gz", "rev": "_R2_001.fastq.gz"},
}


class PonParams:
    """Parameters related to the PON creation workflow."""

    MIN_PON_SAMPLES: int = 6


BIOINFO_TOOL_ENV: Dict[str, str] = {
    BioinfoTools.BEDTOOLS: DockerContainers.ALIGN_QC,
    BioinfoTools.BWA: DockerContainers.ALIGN_QC,
    BioinfoTools.COMPRESS: DockerContainers.ALIGN_QC,
    BioinfoTools.FASTQC: DockerContainers.ALIGN_QC,
    BioinfoTools.SAMTOOLS: DockerContainers.ALIGN_QC,
    BioinfoTools.PICARD: DockerContainers.ALIGN_QC,
    BioinfoTools.MULTIQC: DockerContainers.MULTIQC,
    BioinfoTools.FASTP: DockerContainers.ALIGN_QC,
    BioinfoTools.CSVKIT: DockerContainers.ALIGN_QC,
    BioinfoTools.VEP: DockerContainers.ANNOTATE,
    BioinfoTools.GENMOD: DockerContainers.ANNOTATE,
    BioinfoTools.VCFANNO: DockerContainers.ANNOTATE,
    BioinfoTools.SAMBAMBA: DockerContainers.COVERAGE_QC,
    BioinfoTools.MOSDEPTH: DockerContainers.COVERAGE_QC,
    BioinfoTools.MSISENSORPRO: DockerContainers.MSISENSORPRO,
    BioinfoTools.BCFTOOLS: DockerContainers.PYTHON_3,
    BioinfoTools.TABIX: DockerContainers.PYTHON_3,
    BioinfoTools.BGZIP: DockerContainers.PYTHON_3,
    BioinfoTools.VARDICT: DockerContainers.PYTHON_3,
    BioinfoTools.SVDB: DockerContainers.PYTHON_3,
    BioinfoTools.TIDDIT: DockerContainers.PYTHON_3,
    BioinfoTools.CNVPYTOR: DockerContainers.CNVPYTOR,
    BioinfoTools.MANTA: DockerContainers.PYTHON_27,
    BioinfoTools.CNVKIT: DockerContainers.CNVKIT,
    BioinfoTools.DELLY: DockerContainers.DELLY,
    BioinfoTools.ASCAT: DockerContainers.ASCAT,
    BioinfoTools.VCF2CYTOSURE: DockerContainers.VCF2CYTOSURE,
    BioinfoTools.SOMALIER: DockerContainers.SOMALIER,
    BioinfoTools.CADD: DockerContainers.CADD,
    BioinfoTools.PURECN: DockerContainers.PURECN,
    BioinfoTools.GATK: DockerContainers.GATK,
    BioinfoTools.PYSAM: DockerContainers.PYTHON_3,
}
