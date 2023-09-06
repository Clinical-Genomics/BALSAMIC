"""Balsamic analysis workflow constants."""
from typing import Dict, List

from BALSAMIC.constants.cache import DockerContainers
from BALSAMIC.utils.class_types import StrEnum


class RunMode(StrEnum):
    """Balsamic workflow run mode."""

    CLUSTER: str = "cluster"
    LOCAL: str = "local"


class Gender(StrEnum):
    """Sex options."""

    FEMALE: str = "female"
    MALE: str = "male"


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


class SequencingType(StrEnum):
    """Sequencing carried out."""

    TARGETED: str = "targeted"
    WGS: str = "wgs"


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


class RuleDeliveryMode(StrEnum):
    """Rules to deliver mode."""

    APPEND: str = "append"
    RESET: str = "reset"


RULE_DELIVERY_MODES: List[RuleDeliveryMode] = [mode.value for mode in RuleDeliveryMode]


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
    MOSDEPTH: str = "mosdepth"
    MULTIQC: str = "multiqc"
    PICARD: str = "picard"
    SAMBAMBA: str = "sambamba"
    SAMTOOLS: str = "samtools"
    SOMALIER: str = "somalier"
    SVDB: str = "svdb"
    TABIX: str = "tabix"
    TIDDIT: str = "tiddit"
    VARDICT: str = "vardict"
    VCF2CYTOSURE: str = "vcf2cytosure"
    VCFANNO: str = "vcfanno"
    CADD: str = "cadd"


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
    BioinfoTools.BEDTOOLS.value: DockerContainers.ALIGN_QC.value,
    BioinfoTools.BWA.value: DockerContainers.ALIGN_QC.value,
    BioinfoTools.COMPRESS.value: DockerContainers.ALIGN_QC.value,
    BioinfoTools.FASTQC.value: DockerContainers.ALIGN_QC.value,
    BioinfoTools.SAMTOOLS.value: DockerContainers.ALIGN_QC.value,
    BioinfoTools.PICARD.value: DockerContainers.ALIGN_QC.value,
    BioinfoTools.MULTIQC.value: DockerContainers.ALIGN_QC.value,
    BioinfoTools.FASTP.value: DockerContainers.ALIGN_QC.value,
    BioinfoTools.CSVKIT.value: DockerContainers.ALIGN_QC.value,
    BioinfoTools.VEP.value: DockerContainers.ANNOTATE.value,
    BioinfoTools.GENMOD.value: DockerContainers.ANNOTATE.value,
    BioinfoTools.VCFANNO.value: DockerContainers.ANNOTATE.value,
    BioinfoTools.SAMBAMBA.value: DockerContainers.COVERAGE_QC.value,
    BioinfoTools.MOSDEPTH.value: DockerContainers.COVERAGE_QC.value,
    BioinfoTools.BCFTOOLS.value: DockerContainers.PYTHON_3.value,
    BioinfoTools.TABIX.value: DockerContainers.PYTHON_3.value,
    BioinfoTools.BGZIP.value: DockerContainers.PYTHON_3.value,
    BioinfoTools.GATK.value: DockerContainers.PYTHON_3.value,
    BioinfoTools.VARDICT.value: DockerContainers.PYTHON_3.value,
    BioinfoTools.SVDB.value: DockerContainers.PYTHON_3.value,
    BioinfoTools.TIDDIT.value: DockerContainers.PYTHON_3.value,
    BioinfoTools.CNVPYTOR.value: DockerContainers.CNVPYTOR.value,
    BioinfoTools.MANTA.value: DockerContainers.PYTHON_27.value,
    BioinfoTools.CNVKIT.value: DockerContainers.CNVKIT.value,
    BioinfoTools.DELLY.value: DockerContainers.DELLY.value,
    BioinfoTools.ASCAT.value: DockerContainers.ASCAT.value,
    BioinfoTools.VCF2CYTOSURE.value: DockerContainers.VCF2CYTOSURE.value,
    BioinfoTools.SOMALIER.value: DockerContainers.SOMALIER.value,
    BioinfoTools.CADD.value: DockerContainers.CADD.value,
}
