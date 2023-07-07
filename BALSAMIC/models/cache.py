"""Balsamic reference cache models."""
import logging
from pathlib import Path
from typing import Dict, Optional, List, Any, Union

from pydantic import BaseModel, AnyUrl, DirectoryPath, validator, FilePath

from BALSAMIC.constants.cache import (
    GenomeVersion,
    GRCHVersion,
)
from BALSAMIC.constants.constants import FileType, BwaIndexFileType
from BALSAMIC.utils.exc import BalsamicError

LOG = logging.getLogger(__name__)


class ReferenceUrl(BaseModel):
    """
    Reference model handling URLs and destination paths.

    Attributes:
        url (AnyUrl)              : Address of the reference to download.
        file_type (FileType)      : Reference file type.
        gzip (bool)               : Compression status.
        file_name (str)           : Reference file name after being downloaded.
        dir_name (str)            : Destination directory of the downloaded file.
        file_path (str, optional) : Full reference file path.
        secret (str)              : Database key.

    """

    url: AnyUrl
    file_type: FileType
    gzip: bool
    file_name: str
    dir_name: str
    file_path: Optional[str]
    secret: Optional[str]


class References(BaseModel):
    """
    Reference files model.

    Attributes:
        genome_chrom_size (ReferenceUrl) : Genome chromosome sizes.
        reference_genome (ReferenceUrl)  : Required field for reference genome FASTA file.
        refgene_sql (ReferenceUrl)       : RefSeq's gene SQL format from UCSC.
        refgene_txt (ReferenceUrl)       : RefSeq's gene flat format from UCSC.
    """

    genome_chrom_size: ReferenceUrl
    reference_genome: ReferenceUrl
    refgene_sql: ReferenceUrl
    refgene_txt: ReferenceUrl

    def get_reference_genome_files(self) -> List[str]:
        """Return output reference genome files."""
        return [
            self.reference_genome.file_path,
            self.reference_genome.file_path + "." + FileType.FAI,
            self.reference_genome.file_path.replace(FileType.FASTA, FileType.DICT),
        ] + self.get_reference_genome_bwa_index_files()

    def get_reference_genome_bwa_index_files(self) -> List[str]:
        """Return output BWA reference genome index files."""
        bwa_index_files: List[str] = []
        for bwa_index in BwaIndexFileType:
            bwa_index_files.append(self.reference_genome.file_path + "." + bwa_index)
        return bwa_index_files

    def get_refgene_files(self) -> List[str]:
        """Return RefSeq's gene files from UCSC."""
        return [
            self.refgene_txt.file_path,
            self.get_refgene_flat_file(),
            self.get_refgene_bed_file(),
        ]

    def get_refgene_flat_file(self) -> str:
        """Return RefSeq's gene flat file from UCSC."""
        return self.refgene_txt.file_path.replace(FileType.TXT, FileType.FLAT)

    def get_refgene_bed_file(self) -> str:
        """Return RefSeq's gene BED file from UCSC."""
        return self.get_refgene_flat_file() + "." + FileType.BED


class ReferencesCanFam(References):
    """Canine reference genome files model."""


class ReferencesHg(References):
    """
    Human reference genome files model.

    Attributes:
        access_regions (ReferenceUrl)           : Accessible genome regions.
        ascat_chr_y_loci (ReferenceUrl)         : Chromosome Y loci.
        ascat_gc_correction (ReferenceUrl)      : Genome GC correction bins.
        cadd_snv (ReferenceUrl)                 : CADD SNV annotation file.
        clinvar (ReferenceUrl)                  : ClinVar reference.
        cosmic (ReferenceUrl)                   : COSMIC database's variants as VCF.
        dbsnp (ReferenceUrl)                    : dbSNP VCF file.
        delly_exclusion (ReferenceUrl)          : Genome exclusion regions.
        delly_mappability (ReferenceUrl)        : Genome mappability.
        delly_mappability_findex (ReferenceUrl) : Genome mappability fasta index.
        delly_mappability_gindex (ReferenceUrl) : Genome mappability index.
        gnomad_variant (ReferenceUrl)           : gnomAD variants (non SV) as VCF.
        gnomad_variant_index (ReferenceUrl)     : gnomAD variants VCF index.
        hc_vcf_1kg (ReferenceUrl)               : High confidence 1000 Genome VCF.
        known_indel_1kg (ReferenceUrl)          : 1000 Genome known InDels VCF.
        mills_1kg (ReferenceUrl)                : Mills' high confidence InDels VCF.
        rank_score (ReferenceUrl)               : Rank score model.
        somalier_sites (ReferenceUrl)           : Somalier sites VCF.
        vcf_1kg (ReferenceUrl)                  : 1000 Genome all SNPs.
        wgs_calling_regions (ReferenceUrl)      : WGS calling intervals.
    """

    access_regions: ReferenceUrl
    ascat_chr_y_loci: ReferenceUrl
    ascat_gc_correction: ReferenceUrl
    cadd_snv: ReferenceUrl
    clinvar: ReferenceUrl
    cosmic: ReferenceUrl
    dbsnp: ReferenceUrl
    delly_exclusion: ReferenceUrl
    delly_mappability: ReferenceUrl
    delly_mappability_findex: ReferenceUrl
    delly_mappability_gindex: ReferenceUrl
    gnomad_variant: ReferenceUrl
    gnomad_variant_index: ReferenceUrl
    hc_vcf_1kg: ReferenceUrl
    known_indel_1kg: ReferenceUrl
    mills_1kg: ReferenceUrl
    rank_score: ReferenceUrl
    somalier_sites: ReferenceUrl
    vcf_1kg: ReferenceUrl
    wgs_calling_regions: ReferenceUrl

    def get_cadd_snv_files(self) -> List[str]:
        """Return CADD SNV reference output files."""
        return [self.cadd_snv.file_path, self.cadd_snv.file_path + "." + FileType.TBI]

    def get_delly_files(self) -> List[str]:
        """Return Delly associated output files."""
        return [
            self.delly_exclusion.file_path,
            self.get_delly_exclusion_converted_file(),
            self.delly_mappability.file_path,
            self.delly_mappability_findex.file_path,
            self.delly_mappability_gindex.file_path,
        ]

    def get_delly_exclusion_converted_file_path(self) -> str:
        """Return path to Delly exclusion converted file."""
        return self.delly_exclusion.file_path.replace(
            "." + FileType.TSV, "_converted." + FileType.TSV
        )

    def get_gnomad_files(self) -> List[str]:
        """Return gnomAD associated output files."""
        return [
            self.gnomad_variant.file_path,
            self.gnomad_variant_index.file_path,
        ]

    def get_1k_genome_files(self) -> List[str]:
        """Return 1000 Genome related files."""
        return [
            self.known_indel_1kg.file_path + "." + FileType.GZ,
            self.mills_1kg.file_path + "." + FileType.GZ,
            self.hc_vcf_1kg.file_path + "." + FileType.GZ,
            self.vcf_1kg.file_path + "." + FileType.GZ,
        ]


class AnalysisReferences(BaseModel):
    """
    Reference files pytest t for a general Balsamic analysis.

    Attributes:
        genome_chrom_size (FilePath) : Genome chromosome sizes.
        reference_genome (FilePath)  : Required field for reference genome FASTA file.
        refgene_bed (FilePath)       : RefSeq's gene BED format from UCSC.
        refgene_flat (FilePath)      : RefSeq's gene flat format from UCSC.
        refgene_txt (FilePath)       : RefSeq's gene txt format from UCSC.
    """

    genome_chrom_size: FilePath
    reference_genome: FilePath
    refgene_bed: FilePath
    refgene_flat: FilePath
    refgene_txt: FilePath


class AnalysisReferencesCanFam(AnalysisReferences):
    """Canine reference genome files model."""


class AnalysisReferencesHg(AnalysisReferences):
    """
    Human reference genome files model.

    Attributes:
        access_regions (FilePath)            : Accessible genome regions.
        ascat_chr_y_loci (FilePath)          : Chromosome Y loci.
        ascat_gc_correction (FilePath)       : Genome GC correction bins.
        cadd_snv (FilePath)                  : CADD SNV annotation file.
        clinvar (FilePath)                   : ClinVar reference.
        cosmic (FilePath)                    : COSMIC database's variants as VCF.
        dbsnp (FilePath)                     : dbSNP VCF file.
        delly_exclusion (FilePath)           : Genome exclusion regions.
        delly_exclusion_converted (FilePath) : Genome exclusion regions without "chr" field.
        delly_mappability (FilePath)         : Genome mappability.
        gnomad_variant (FilePath)            : gnomAD variants (non SV) as VCF.
        hc_vcf_1kg (FilePath)                : High confidence 1000 Genome VCF.
        known_indel_1kg (FilePath)           : 1000 Genome known InDels VCF.
        mills_1kg (FilePath)                 : Mills' high confidence InDels VCF.
        rank_score (FilePath)                : Rank score model.
        somalier_sites (FilePath)            : Somalier sites VCF.
        vcf_1kg (FilePath)                   : 1000 Genome all SNPs.
        vep_dir (DirectoryPath)              : VEP annotations output directory.
        wgs_calling_regions (FilePath)       : WGS calling intervals.
    """

    access_regions: FilePath
    ascat_chr_y_loci: FilePath
    ascat_gc_correction: FilePath
    cadd_snv: FilePath
    clinvar: FilePath
    cosmic: FilePath
    dbsnp: FilePath
    delly_exclusion: FilePath
    delly_exclusion_converted: FilePath
    delly_mappability: FilePath
    gnomad_variant: FilePath
    hc_vcf_1kg: FilePath
    known_indel_1kg: FilePath
    mills_1kg: FilePath
    rank_score: FilePath
    somalier_sites: FilePath
    vcf_1kg: FilePath
    vep_dir: DirectoryPath
    wgs_calling_regions: FilePath


class CacheAnalysis(BaseModel):
    """
    Reference analysis configuration model.


    Attributes:
        case_id (str) : Reference case identifier.
    """

    case_id: str


class CacheConfig(BaseModel):
    """
    Reference build configuration model.

    Attributes:
        analysis (CacheAnalysis)                           : Reference analysis model.
        references_dir (DirectoryPath)                     : Output directory for the downloaded reference.
        containers_dir (DirectoryPath)                     : Output directory for the downloaded singularity containers.
        genome_dir (DirectoryPath)                         : Genome references output directory.
        variants_dir (DirectoryPath)                       : Variant references output directory for the.
        vep_dir (Path)                                     : VEP annotations output directory.
        genome_version (GenomeVersion)                     : Genome version associated with the balsamic cache.
        cosmic_key (str, optional)                         : COSMIC database key.
        bioinfo_tools (dict)                               : Dictionary of bioinformatics software and containers.
        containers (Dict[str, str])                        : Dictionary linking container names and dockerhub images.
        references (Union[ReferencesHg, ReferencesCanFam]) : Reference files model.
        references_date (str)                              : Reference access date.
    """

    analysis: CacheAnalysis
    references_dir: DirectoryPath
    containers_dir: DirectoryPath
    genome_dir: DirectoryPath
    variants_dir: DirectoryPath
    vep_dir: Path
    genome_version: GenomeVersion
    cosmic_key: Optional[str]
    bioinfo_tools: dict
    containers: Dict[str, str]
    references: Union[ReferencesHg, ReferencesCanFam]
    references_date: str

    @validator("references")
    def validate_references(
        cls, references: References, values: Dict[str, Any]
    ) -> References:
        """Validate the reference output paths."""
        for model in references:
            reference_key: str
            reference: ReferenceUrl
            reference_key, reference = model[0], model[1]
            reference.file_path = (
                Path(
                    values.get("references_dir"),
                    reference.dir_name,
                    reference.file_name,
                ).as_posix()
                if reference
                else None
            )
            reference.secret = (
                values.get("cosmic_key") if "cosmic" in reference_key else None
            )
        return references

    def get_grch_version(self) -> Optional[GRCHVersion]:
        """Return GRCH format version of the genome version."""
        version: Dict[GenomeVersion, GRCHVersion] = {
            GenomeVersion.HG19: GRCHVersion.GRCH37,
            GenomeVersion.HG38: GRCHVersion.GRCH38,
        }
        return version.get(self.genome_version)

    def get_reference_paths(self) -> List[str]:
        """Return a list of reference paths."""
        return [
            Path(reference[1].file_path).as_posix() for reference in self.references
        ]

    def get_reference_by_path(self, reference_path: str) -> ReferenceUrl:
        """Return a reference given its full path."""
        for model in self.references:
            reference: ReferenceUrl = model[1]
            if reference.file_path == reference_path:
                return reference
        LOG.error(f"No reference with the provided reference path {reference_path}")
        raise BalsamicError()

    def get_reference_paths_by_file_type_and_compression(
        self, file_type: FileType, compression: bool
    ) -> List[str]:
        """Return a list of reference paths given a file type and a compression status."""
        file_type_references: List[str] = self.get_reference_paths_by_file_type(
            file_type=file_type
        )
        compression_references: List[str] = self.get_reference_paths_by_compression(
            compression=compression
        )
        return list(set(file_type_references).intersection(compression_references))

    def get_reference_paths_by_file_type(self, file_type: FileType) -> List[str]:
        """Return a list of reference paths given a file type."""
        return [
            reference[1].file_path
            for reference in self.references
            if reference[1].file_type == file_type
        ]

    def get_reference_paths_by_compression(self, compression: bool) -> List[str]:
        """Return a list of reference paths given a compression status."""
        return [
            reference[1].file_path
            for reference in self.references
            if reference[1].gzip == compression
        ]

    def get_compressed_indexed_vcfs(self) -> List[str]:
        """Return an output list of compressed and indexed VCFs."""
        return [
            vcf + ".gz.tbi"
            for vcf in self.get_reference_paths_by_file_type_and_compression(
                file_type=FileType.VCF, compression=True
            )
        ]

    def get_container_output_paths(self) -> List[str]:
        """Return a complete list of output singularity images."""
        return [
            Path(self.containers_dir, image + "." + FileType.SIF).as_posix()
            for image in self.containers.keys()
        ]

    def get_reference_output_paths(self) -> List[str]:
        """Return a complete list of output reference paths."""
        reference_paths: List[str] = [
            self.references.genome_chrom_size.file_path,
            *self.references.get_reference_genome_files(),
            *self.references.get_refgene_files(),
        ]
        if self.genome_version == GenomeVersion.CanFam3:
            return reference_paths

        reference_paths += [
            self.references.access_regions.file_path,
            self.references.ascat_chr_y_loci.file_path,
            self.references.ascat_gc_correction.file_path,
            self.references.cadd_snv.file_path,
            self.references.clinvar.file_path + "." + FileType.GZ,
            self.references.cosmic.file_path + "." + FileType.GZ,
            self.references.dbsnp.file_path + "." + FileType.GZ,
            self.references.rank_score.file_path,
            self.references.somalier_sites.file_path + "." + FileType.GZ,
            self.references.wgs_calling_regions.file_path,
            *self.get_compressed_indexed_vcfs(),
            *self.references.get_1k_genome_files(),
            *self.references.get_cadd_snv_files(),
            *self.references.get_delly_files(),
            *self.references.get_gnomad_files(),
            self.vep_dir.as_posix(),
        ]
        return reference_paths

    def get_analysis_references(
        self,
    ) -> Union[AnalysisReferencesHg, AnalysisReferencesCanFam]:
        """Return reference output model for Balsamic analyses."""
        if self.genome_version == GenomeVersion.CanFam3:
            return AnalysisReferencesCanFam(
                genome_chrom_size=self.references.genome_chrom_size.file_path,
                reference_genome=self.references.reference_genome.file_path,
                refgene_bed=self.references.get_refgene_bed_file(),
                refgene_flat=self.references.get_refgene_flat_file(),
                refgene_txt=self.references.refgene_txt.file_path,
            )

        return AnalysisReferencesHg(
            access_regions=self.references.access_regions.file_path,
            ascat_chr_y_loci=self.references.ascat_chr_y_loci.file_path,
            ascat_gc_correction=self.references.ascat_gc_correction.file_path,
            cadd_snv=self.references.cadd_snv.file_path,
            clinvar=self.references.clinvar.file_path + "." + FileType.GZ,
            cosmic=self.references.cosmic.file_path + "." + FileType.GZ,
            dbsnp=self.references.dbsnp.file_path + "." + FileType.GZ,
            delly_exclusion=self.references.delly_exclusion.file_path,
            delly_exclusion_converted=self.references.get_delly_exclusion_converted_file(),
            delly_mappability=self.references.delly_mappability.file_path,
            genome_chrom_size=self.references.genome_chrom_size.file_path,
            gnomad_variant=self.references.gnomad_variant.file_path,
            hc_vcf_1kg=self.references.hc_vcf_1kg.file_path + "." + FileType.GZ,
            known_indel_1kg=self.references.known_indel_1kg.file_path
            + "."
            + FileType.GZ,
            mills_1kg=self.references.mills_1kg.file_path + "." + FileType.GZ,
            rank_score=self.references.rank_score.file_path,
            reference_genome=self.references.reference_genome.file_path,
            refgene_bed=self.references.get_refgene_bed_file(),
            refgene_flat=self.references.get_refgene_flat_file(),
            refgene_txt=self.references.refgene_txt.file_path,
            somalier_sites=self.references.somalier_sites.file_path + "." + FileType.GZ,
            vcf_1kg=self.references.vcf_1kg.file_path + "." + FileType.GZ,
            vep_dir=self.vep_dir.as_posix(),
            wgs_calling_regions=self.references.wgs_calling_regions.file_path,
        )
