"""Balsamic cache models."""
import hashlib
import logging
from pathlib import Path
from typing import Dict, Optional, List, Any

from pydantic import BaseModel, AnyUrl, DirectoryPath, validator, FilePath

from BALSAMIC.constants.analysis import GenomeVersion
from BALSAMIC.constants.cache import FileType, BwaIndexFileType

LOG = logging.getLogger(__name__)


class ReferenceUrlModel(BaseModel):
    """
    Reference model handling URLs.

    Attributes:
        url       : address of the reference to download
        file_type : reference file type
        gzip      : compression status
        file_name : reference file name after being downloaded
        dir_name  : destination directory of the downloaded file
        file_path : full reference file path


    """

    url: AnyUrl
    file_type: FileType
    gzip: bool
    file_name: str
    dir_name: str
    file_path: Optional[str]

    @property
    def write_md5(self) -> None:
        """Write the md5 checksum for the first 4 kB of a file."""
        if not Path(file_path).is_file():
            LOG.error(f"File {self.file_path} does not exist")
            raise FileNotFoundError
        with open(self.file_path, "rb") as fh:
            for chunk in iter(lambda: fh.read(4096), b""):
                hashlib.md5().update(chunk)
        with open(file_path + ".md5", "w") as fh:
            fh.write("{} {}\n".format(file_path, hashlib.md5().hexdigest()))


class ReferencesModel(BaseModel):
    """
    Reference files model.

    Attributes:
        access_regions           : accessible genome regions
        ascat_chryloci           : chromosome Y loci
        ascat_gccorrection       : genome GC correction bins
        clinvar                  : ClinVar reference
        cosmicdb                 : COSMIC database's variants as VCF
        dbsnp                    : dbSNP VCF file
        delly_exclusion          : genome exclusion regions
        delly_mappability        : genome mappability
        delly_mappability_findex : genome mappability fasta index
        delly_mappability_gindex : genome mappability index
        genome_chrom_size        : genome chromosome sizes
        gnomad_variant           : gnomAD variants (non SV) as VCF
        gnomad_variant_index     : gnomAD variants VCF index
        hc_vcf_1kg               : high confidence 1000 Genome VCF
        known_indel_1kg          : 1000 Genome known InDels VCF
        mills_1kg                : Mills' high confidence InDels VCF
        rankscore                : rank score model
        reference_genome         : required field for reference genome FASTA file
        refgene_sql              : RefSeq's gene SQL format from UCSC
        refgene_txt              : RefSeq's gene flat format from UCSC
        somalier_sites           : somalier sites VCF
        vcf_1kg                  : 1000 Genome all SNPs
        wgs_calling              : WGS calling intervals
    """

    access_regions: Optional[ReferenceUrlModel]
    ascat_chryloci: Optional[ReferenceUrlModel]
    ascat_gccorrection: Optional[ReferenceUrlModel]
    clinvar: Optional[ReferenceUrlModel]
    cosmicdb: Optional[ReferenceUrlModel]
    dbsnp: Optional[ReferenceUrlModel]
    delly_exclusion: Optional[ReferenceUrlModel]
    delly_mappability: Optional[ReferenceUrlModel]
    delly_mappability_findex: Optional[ReferenceUrlModel]
    delly_mappability_gindex: Optional[ReferenceUrlModel]
    genome_chrom_size: Optional[ReferenceUrlModel]
    gnomad_variant: Optional[ReferenceUrlModel]
    gnomad_variant_index: Optional[ReferenceUrlModel]
    hc_vcf_1kg: Optional[ReferenceUrlModel]
    known_indel_1kg: Optional[ReferenceUrlModel]
    mills_1kg: Optional[ReferenceUrlModel]
    rankscore: Optional[ReferenceUrlModel]
    reference_genome: ReferenceUrlModel
    refgene_sql: Optional[ReferenceUrlModel]
    refgene_txt: Optional[ReferenceUrlModel]
    somalier_sites: Optional[ReferenceUrlModel]
    vcf_1kg: Optional[ReferenceUrlModel]
    wgs_calling: Optional[ReferenceUrlModel]

    def get_bwa_index_files(self) -> List[str]:
        """Return output BWA genome index files."""
        bwa_index_files: List[str] = []
        for bwa_index in BwaIndexFileType:
            bwa_index_files.append(self.reference_genome.file_path + "." + bwa_index)
        return bwa_index_files

    def get_delly_files(self) -> List[str]:
        """Return delly associated output files."""
        return [
            self.delly_exclusion.file_path,
            self.delly_exclusion.file_path.replace(
                "." + FileType.TSV, "_converted." + FileType.TSV
            ),
            self.delly_mappability.file_path,
            self.delly_mappability_findex.file_path,
            self.delly_mappability_gindex.file_path,
        ]

    def get_gnomad_files(self) -> List[str]:
        """Return gnomad associated output files."""
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

    def get_reference_genome_files(self) -> List[str]:
        """Return output reference genome files."""
        return [
            self.reference_genome.file_path,
            self.reference_genome.file_path + "." + FileType.FAI,
            self.reference_genome.file_path.replace(FileType.FASTA, FileType.DICT),
        ]

    def get_refgene_files(self) -> List[str]:
        """Return RefSeq's gene files from UCSC."""
        return [
            self.refgene_txt.file_path,
            self.refgene_txt.file_path.replace(FileType.TXT, FileType.FLAT),
            self.refgene_txt.file_path.replace(FileType.TXT, FileType.FLAT)
            + "."
            + FileType.BED,
        ]


class CacheConfigModel(BaseModel):
    """
    Reference build configuration model.

    Attributes:
        references_dir  : output directory for the downloaded reference
        containers_dir  : output directory for the downloaded singularity containers
        genome_version  : genome version associated with the balsamic cache
        cosmic_key      : COSMIC database key
        bioinfo_tools   : dictionary of bioinformatics software and their associated containers
        containers      : dictionary linking the container names and their dockerhub image paths
        references      : reference files model
        references_date : reference access date
    """

    references_dir: DirectoryPath
    containers_dir: DirectoryPath
    genome_version: GenomeVersion
    cosmic_key: Optional[str]
    bioinfo_tools: dict
    containers: Dict[str, str]
    references: ReferencesModel
    references_date: str

    @validator("references")
    def validate_reference_output_path(
        cls, references: ReferencesModel, values: Dict[str, Any]
    ) -> ReferencesModel:
        for model in references:
            reference: ReferenceUrlModel = model[1]
            if reference:
                reference.file_path = Path(
                    values.get("references_dir"),
                    reference.dir_name,
                    reference.file_name,
                ).as_posix()
        return references

    def get_reference_paths(self) -> List[str]:
        """Return a list of output reference paths."""
        reference_paths: List[str] = [
            self.references.access_regions.file_path,
            self.references.ascat_chryloci.file_path,
            self.references.ascat_gccorrection.file_path,
            self.references.clinvar.file_path + "." + FileType.GZ,
            self.references.cosmicdb.file_path + "." + FileType.GZ,
            self.references.dbsnp.file_path + "." + FileType.GZ,
            self.references.genome_chrom_size.file_path,
            self.references.rankscore.file_path,
            self.references.somalier_sites.file_path + "." + FileType.GZ,
            self.references.wgs_calling.file_path,
        ]
        reference_paths.extend(
            self.references.get_bwa_index_files()
            + self.references.get_delly_files()
            + self.references.get_gnomad_files()
            + self.references.get_1k_genome_files()
            + self.references.get_reference_genome_files()
            + self.references.get_refgene_files()
            + [
                vcf + ".gz.tbi"
                for vcf in self.get_reference_paths_by_file_type(
                    file_type=FileType.VCF, compression=True
                )
            ]
        )
        return reference_paths

    def get_reference_paths_by_file_type(
        self, file_type: FileType, compression: bool
    ) -> List[str]:
        """Return a list of reference paths given a file type and a compression status."""
        reference_paths: List[str] = []
        for model in self.references:
            reference: ReferenceUrlModel = model[1]
            if reference.file_type == file_type and reference.gzip == compression:
                reference_paths.append(reference.file_path)
        return reference_paths