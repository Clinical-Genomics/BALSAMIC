"""Balsamic reference cache models."""
import hashlib
import logging
from pathlib import Path
from typing import Optional

from pydantic import BaseModel, validator, AnyUrl

from BALSAMIC.constants.cache import GenomeVersion, FileType

LOG = logging.getLogger(__name__)


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
    genome_version: Optional[str]
    file_name: Optional[str]
    dir_name: Optional[str]
    secret: Optional[str]

    @validator("file_type")
    def check_file_type(cls, value) -> str:
        """Validate file format according to constants"""
        assert value in set(FileType), f"{value} not a valid reference file format."
        return value

    @validator("genome_version")
    def check_genome_ver(cls, value) -> str:
        """Validate genome version according constants"""
        assert value in set(GenomeVersion), f"{value} not a valid genome version."
        return value

    @property
    def get_output_file(self):
        """return output file full path"""
        output_file_path = Path(self.dir_name, self.file_name).as_posix()
        return output_file_path

    @property
    def write_md5(self):
        """calculate md5 for first 4kb of file and write to file_name.md5"""
        hash_md5 = hashlib.md5()
        output_file = Path(self.dir_name, self.file_name)
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
        wgs_calling_regions: ReferenceUrlsModel. Optional field for wgs calling intervals
        genome_chrom_size: ReferenceUrlsModel. Optional field for geneome's chromosome sizes
        gnomad_variant: ReferenceUrlsModel. Optional gnomad variants (non SV) as vcf
        cosmic: ReferenceUrlsModel. Optional COSMIC database's variants as vcf
        refgene_txt: ReferenceUrlsModel. Optional refseq's gene flat format from UCSC
        refgene_sql: ReferenceUrlsModel. Optional refseq's gene sql format from UCSC
        rank_score: ReferenceUrlsModel. Optional rankscore model
        access_regions: ReferenceUrlsModel. Optional field for accessible genome regions
        delly_exclusion: ReferenceUrlsModel. Optional field for genome exclusion regions
        delly_mappability: ReferenceUrlsModel. Optional field for genome mappability
        ascat_gc_correction: ReferenceUrlsModel. Optional field for genome gc correction bins
        ascat_chr_y_loci: ReferenceUrlsModel. Optional field for chromosome Y loci
        clinvar: ReferenceUrlsModel. Optional field for clinvar reference
        somalier_sites: ReferenceUrlsModel. Optional field for somalier sites vcf
        cadd_snv: ReferenceUrlsModel. Optional field for CADD SNV reference
    """

    basedir: str = ""
    reference_genome: ReferenceUrlsModel
    dbsnp: Optional[ReferenceUrlsModel]
    hc_vcf_1kg: Optional[ReferenceUrlsModel]
    mills_1kg: Optional[ReferenceUrlsModel]
    known_indel_1kg: Optional[ReferenceUrlsModel]
    vcf_1kg: Optional[ReferenceUrlsModel]
    wgs_calling_regions: Optional[ReferenceUrlsModel]
    genome_chrom_size: Optional[ReferenceUrlsModel]
    gnomad_variant: Optional[ReferenceUrlsModel]
    gnomad_variant_index: Optional[ReferenceUrlsModel]
    cosmic: Optional[ReferenceUrlsModel]
    refgene_txt: Optional[ReferenceUrlsModel]
    refgene_sql: Optional[ReferenceUrlsModel]
    rank_score: Optional[ReferenceUrlsModel]
    access_regions: Optional[ReferenceUrlsModel]
    delly_exclusion: Optional[ReferenceUrlsModel]
    delly_mappability: Optional[ReferenceUrlsModel]
    delly_mappability_gindex: Optional[ReferenceUrlsModel]
    delly_mappability_findex: Optional[ReferenceUrlsModel]
    ascat_gc_correction: Optional[ReferenceUrlsModel]
    ascat_chr_y_loci: Optional[ReferenceUrlsModel]
    clinvar: Optional[ReferenceUrlsModel]
    somalier_sites: Optional[ReferenceUrlsModel]
    cadd_snv: Optional[ReferenceUrlsModel]
    cadd_snv_index: Optional[ReferenceUrlsModel]

    @validator("*", pre=True)
    def validate_path(cls, value, values, **kwargs):
        """validate and append path in ReferenceUrlsModel fields with basedir"""
        if isinstance(value, str):
            output_value = value
        else:
            if "dir_name" in value:
                value["dir_name"] = Path(
                    values["basedir"], value["dir_name"]
                ).as_posix()
                output_value = ReferenceUrlsModel.parse_obj(value)
            else:
                output_value = value

        return output_value
