"""Balsamic cache specific constants."""
from typing import Dict

from BALSAMIC.constants.constants import FileType
from BALSAMIC.utils.class_types import StrEnum

DOCKER_URL: str = "docker://clinicalgenomics/balsamic"
VEP_PLUGINS: str = "all"


class GenomeVersion(StrEnum):
    """Reference genome versions."""

    HG19: str = "hg19"
    HG38: str = "hg38"
    CanFam3: str = "canfam3"


class GRCHVersion(StrEnum):
    """Genome Reference Consortium Human Reference versions."""

    GRCH37: str = "GRCh37"
    GRCH38: str = "GRCh38"


class Species(StrEnum):
    """A class representing different species."""

    HOMO_SAPIENS: str = "homo_sapiens_merged"


class ContainerVersion(StrEnum):
    """Balsamic container versions."""

    DEVELOP: str = "develop"
    RELEASE: str = "release"


class DockerContainers(StrEnum):
    """Docker containers names."""

    ALIGN_QC: str = "align_qc"
    ANNOTATE: str = "annotate"
    ASCAT: str = "ascatNgs"
    CADD: str = "cadd"
    CNVKIT: str = "varcall_cnvkit"
    CNVPYTOR: str = "cnvpytor"
    COVERAGE_QC: str = "coverage_qc"
    DELLY: str = "delly"
    PYTHON_3: str = "varcall_py3"
    PYTHON_27: str = "varcall_py27"
    SOMALIER: str = "somalier"
    VCF2CYTOSURE: str = "vcf2cytosure"


REFERENCE_FILES: Dict[GenomeVersion, Dict[str, dict]] = {
    GenomeVersion.HG19: {
        "reference_genome": {
            "url": "gs://gatk-legacy-bundles/b37/human_g1k_v37.fasta.gz",
            "file_type": FileType.FASTA,
            "gzip": True,
            "file_name": "human_g1k_v37.fasta",
            "dir_name": "genome",
        },
        "dbsnp": {
            "url": "gs://gatk-legacy-bundles/b37/dbsnp_138.b37.vcf.gz",
            "file_type": FileType.VCF,
            "gzip": True,
            "file_name": "dbsnp_grch37_b138.vcf",
            "dir_name": "variants",
        },
        "hc_vcf_1kg": {
            "url": "gs://gatk-legacy-bundles/b37/1000G_phase1.snps.high_confidence.b37.vcf.gz",
            "file_type": FileType.VCF,
            "gzip": True,
            "file_name": "1kg_phase1_snps_high_confidence_b37.vcf",
            "dir_name": "variants",
        },
        "mills_1kg": {
            "url": "gs://gatk-legacy-bundles/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz",
            "file_type": FileType.VCF,
            "gzip": True,
            "file_name": "mills_1kg_index.vcf",
            "dir_name": "variants",
        },
        "known_indel_1kg": {
            "url": "gs://gatk-legacy-bundles/b37/1000G_phase1.indels.b37.vcf.gz",
            "file_type": FileType.VCF,
            "gzip": True,
            "file_name": "1kg_known_indels_b37.vcf",
            "dir_name": "variants",
        },
        "vcf_1kg": {
            "url": "gs://genomics-public-data/ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz",
            "file_type": FileType.VCF,
            "gzip": True,
            "file_name": "1k_genome_wgs_p1_v3_all_sites.vcf",
            "dir_name": "variants",
        },
        "gnomad_variant": {
            "url": "gs://gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.vcf.bgz",
            "file_type": FileType.VCF,
            "gzip": False,
            "file_name": "gnomad.genomes.r2.1.1.sites.vcf.bgz",
            "dir_name": "variants",
        },
        "gnomad_variant_index": {
            "url": "gs://gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.vcf.bgz.tbi",
            "file_type": FileType.VCF,
            "gzip": False,
            "file_name": "gnomad.genomes.r2.1.1.sites.vcf.bgz.tbi",
            "dir_name": "variants",
        },
        "cosmic": {
            "url": "https://cancer.sanger.ac.uk/cosmic/file_download/GRCh37/cosmic/v97/VCF/CosmicCodingMuts.vcf.gz",
            "file_type": FileType.VCF,
            "gzip": True,
            "file_name": "cosmic_coding_muts_v97.vcf",
            "dir_name": "variants",
        },
        "wgs_calling_regions": {
            "url": "gs://gatk-legacy-bundles/b37/wgs_calling_regions.v1.interval_list",
            "file_type": FileType.TXT,
            "gzip": False,
            "file_name": "wgs_calling_regions.v1",
            "dir_name": "genome",
        },
        "genome_chrom_size": {
            "url": "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes",
            "file_type": FileType.TXT,
            "gzip": False,
            "file_name": "hg19.chrom.sizes",
            "dir_name": "genome",
        },
        "refgene_txt": {
            "url": "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz",
            "file_type": FileType.TXT,
            "gzip": True,
            "file_name": "refGene.txt",
            "dir_name": "genome",
        },
        "refgene_sql": {
            "url": "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.sql",
            "file_type": FileType.TXT,
            "gzip": False,
            "file_name": "refGene.sql",
            "dir_name": "genome",
        },
        "rank_score": {
            "url": "https://raw.githubusercontent.com/Clinical-Genomics/reference-files/master/cancer/rank_model/cancer_rank_model_-v0.1-.ini",
            "file_type": FileType.TXT,
            "gzip": False,
            "file_name": "cancer_rank_model_-v0.1-.ini",
            "dir_name": "genome",
        },
        "access_regions": {
            "url": "https://raw.githubusercontent.com/etal/cnvkit/master/data/access-5k-mappable.hg19.bed",
            "file_type": FileType.TXT,
            "gzip": False,
            "file_name": "access_5kb_hg19.txt",
            "dir_name": "genome",
        },
        "delly_exclusion": {
            "url": "https://raw.githubusercontent.com/dellytools/delly/master/excludeTemplates/human.hg19.excl.tsv",
            "file_type": FileType.TSV,
            "gzip": False,
            "file_name": "delly_exclusion.tsv",
            "dir_name": "genome",
        },
        "delly_mappability": {
            "url": "https://raw.githubusercontent.com/Clinical-Genomics/reference-files/86aab2d10c5ffc009bc8c68ad077ab7283d8fe06/cancer/references/GRCh37.delly.blacklist.gz",
            "file_type": FileType.TXT,
            "gzip": False,
            "file_name": "delly_mappability.gz",
            "dir_name": "genome",
        },
        "delly_mappability_gindex": {
            "url": "https://raw.githubusercontent.com/Clinical-Genomics/reference-files/86aab2d10c5ffc009bc8c68ad077ab7283d8fe06/cancer/references/GRCh37.delly.blacklist.gz.gzi",
            "file_type": FileType.TXT,
            "gzip": False,
            "file_name": "delly_mappability.gz.gzi",
            "dir_name": "genome",
        },
        "delly_mappability_findex": {
            "url": "https://raw.githubusercontent.com/Clinical-Genomics/reference-files/86aab2d10c5ffc009bc8c68ad077ab7283d8fe06/cancer/references/GRCh37.delly.blacklist.gz.fai",
            "file_type": FileType.TXT,
            "gzip": False,
            "file_name": "delly_mappability.gz.fai",
            "dir_name": "genome",
        },
        "ascat_gc_correction": {
            "url": "https://raw.githubusercontent.com/Clinical-Genomics/reference-files/12a6c760fd542c02de2cda286b6245e46f4b6a97/cancer/references/GRCh37_SnpGcCorrections.tsv.gz",
            "file_type": FileType.TSV,
            "gzip": True,
            "file_name": "GRCh37_SnpGcCorrections.tsv",
            "dir_name": "genome",
        },
        "ascat_chr_y_loci": {
            "url": "https://raw.githubusercontent.com/Clinical-Genomics/reference-files/12a6c760fd542c02de2cda286b6245e46f4b6a97/cancer/references/GRCh37_d5_Y.loci",
            "file_type": FileType.TXT,
            "gzip": False,
            "file_name": "GRCh37_Y.loci",
            "dir_name": "genome",
        },
        "clinvar": {
            "url": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz",
            "file_type": FileType.VCF,
            "gzip": True,
            "file_name": "clinvar.vcf",
            "dir_name": "variants",
        },
        "somalier_sites": {
            "url": "https://raw.githubusercontent.com/Clinical-Genomics/reference-files/87f22d3f458569afbcb4d7f1588468d21d1751fb/cancer/references/GRCh37.somalier.sites.vcf.gz",
            "file_type": FileType.VCF,
            "gzip": True,
            "file_name": "GRCh37.somalier.sites.vcf",
            "dir_name": "variants",
        },
        "cadd_snv": {
            "url": "https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh37/whole_genome_SNVs.tsv.gz",
            "file_type": FileType.TSV,
            "gzip": False,
            "file_name": "hg19.cadd_snv.tsv.gz",
            "dir_name": "variants",
        },
    },
    GenomeVersion.HG38: {
        "reference_genome": {
            "url": "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta",
            "file_type": FileType.FASTA,
            "gzip": False,
            "file_name": "Homo_sapiens_assembly38.fasta",
            "dir_name": "genome",
        },
        "dbsnp": {
            "url": "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf",
            "file_type": FileType.VCF,
            "gzip": False,
            "file_name": "Homo_sapiens_assembly38.dbsnp138.vcf",
            "dir_name": "variants",
        },
        "hc_vcf_1kg": {
            "url": "gs://genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
            "file_type": FileType.VCF,
            "gzip": True,
            "file_name": "1000G_phase1.snps.high_confidence.hg38.vcf",
            "dir_name": "variants",
        },
        "mills_1kg": {
            "url": "gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
            "file_type": FileType.VCF,
            "gzip": True,
            "file_name": "Mills_and_1000G_gold_standard.indels.hg38.vcf",
            "dir_name": "variants",
        },
        "known_indel_1kg": {
            "url": "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz",
            "file_type": FileType.VCF,
            "gzip": True,
            "file_name": "Homo_sapiens_assembly38.known_indels.vcf",
            "dir_name": "variants",
        },
        "vcf_1kg": {
            "url": "gs://genomics-public-data/resources/broad/hg38/v0/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf",
            "file_type": FileType.VCF,
            "gzip": False,
            "file_name": "1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf",
            "dir_name": "variants",
        },
        "gnomad_variant": {
            "url": "gs://gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz",
            "file_type": FileType.VCF,
            "gzip": False,
            "file_name": "gnomad.genomes.r2.1.1.sites.vcf.bgz",
            "dir_name": "variants",
        },
        "gnomad_variant_index": {
            "url": "gs://gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz.tbi",
            "file_type": FileType.VCF,
            "gzip": False,
            "file_name": "gnomad.genomes.r2.1.1.sites.vcf.bgz.tbi",
            "dir_name": "variants",
        },
        "cosmic": {
            "url": "https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v97/VCF/CosmicCodingMuts.vcf.gz",
            "file_type": FileType.VCF,
            "gzip": True,
            "file_name": "cosmic_coding_muts_v97.vcf",
            "dir_name": "variants",
        },
        "wgs_calling_regions": {
            "url": "gs://genomics-public-data/resources/broad/hg38/v0/wgs_calling_regions.hg38.interval_list",
            "file_type": FileType.TXT,
            "gzip": False,
            "file_name": "wgs_calling_regions.v1",
            "dir_name": "genome",
        },
        "genome_chrom_size": {
            "url": "https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes",
            "file_type": FileType.TXT,
            "gzip": False,
            "file_name": "hg38.chrom.sizes",
            "dir_name": "genome",
        },
        "refgene_txt": {
            "url": "https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz",
            "file_type": FileType.TXT,
            "gzip": True,
            "file_name": "refGene.txt",
            "dir_name": "genome",
        },
        "refgene_sql": {
            "url": "https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.sql",
            "file_type": FileType.TXT,
            "gzip": False,
            "file_name": "refGene.sql",
            "dir_name": "genome",
        },
        "rank_score": {
            "url": "https://raw.githubusercontent.com/Clinical-Genomics/reference-files/master/cancer/rank_model/cancer_rank_model_-v0.1-.ini",
            "file_type": FileType.TXT,
            "gzip": False,
            "file_name": "cancer_rank_model_-v0.1-.ini",
            "dir_name": "genome",
        },
        "access_regions": {
            "url": "https://raw.githubusercontent.com/etal/cnvkit/master/data/access-5k-mappable.hg19.bed",
            "file_type": FileType.TXT,
            "gzip": False,
            "file_name": "access_5kb_hg38.txt",
            "dir_name": "genome",
        },
        "delly_exclusion": {
            "url": "https://raw.githubusercontent.com/dellytools/delly/master/excludeTemplates/human.hg38.excl.tsv",
            "file_type": FileType.TSV,
            "gzip": False,
            "file_name": "delly_exclusion.tsv",
            "dir_name": "genome",
        },
        "delly_mappability": {
            "url": "https://raw.githubusercontent.com/Clinical-Genomics/reference-files/ea051b864d18945980f0ded6b16a5d192bd736a5/cancer/references/GRCh38.delly.blacklist.gz",
            "file_type": FileType.TXT,
            "gzip": False,
            "file_name": "delly_mappability.gz",
            "dir_name": "genome",
        },
        "delly_mappability_gindex": {
            "url": "https://raw.githubusercontent.com/Clinical-Genomics/reference-files/ea051b864d18945980f0ded6b16a5d192bd736a5/cancer/references/GRCh38.delly.blacklist.gz.gzi",
            "file_type": FileType.TXT,
            "gzip": False,
            "file_name": "delly_mappability.gz.gzi",
            "dir_name": "genome",
        },
        "delly_mappability_findex": {
            "url": "https://raw.githubusercontent.com/Clinical-Genomics/reference-files/ea051b864d18945980f0ded6b16a5d192bd736a5/cancer/references/GRCh38.delly.blacklist.gz.fai",
            "file_type": FileType.TXT,
            "gzip": False,
            "file_name": "delly_mappability.gz.fai",
            "dir_name": "genome",
        },
        "ascat_gc_correction": {
            "url": "https://raw.githubusercontent.com/Clinical-Genomics/reference-files/35465e2644f76f2d59427a9b379d34ecea71f259/cancer/references/hg38_SnpGcCorrections.tsv.gz",
            "file_type": FileType.TSV,
            "gzip": True,
            "file_name": "hg38_SnpGcCorrections.tsv",
            "dir_name": "genome",
        },
        "ascat_chr_y_loci": {
            "url": "https://raw.githubusercontent.com/Clinical-Genomics/reference-files/ea051b864d18945980f0ded6b16a5d192bd736a5/cancer/references/hg38_Y.loci",
            "file_type": FileType.TXT,
            "gzip": False,
            "file_name": "hg38_Y.loci",
            "dir_name": "genome",
        },
        "clinvar": {
            "url": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz",
            "file_type": FileType.VCF,
            "gzip": True,
            "file_name": "clinvar.vcf",
            "dir_name": "variants",
        },
        "somalier_sites": {
            "url": "https://raw.githubusercontent.com/Clinical-Genomics/reference-files/87f22d3f458569afbcb4d7f1588468d21d1751fb/cancer/references/hg38.somalier.sites.vcf.gz",
            "file_type": FileType.VCF,
            "gzip": True,
            "file_name": "hg38.somalier.sites.vcf",
            "dir_name": "variants",
        },
        "cadd_snv": {
            "url": "https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz",
            "file_type": FileType.TSV,
            "gzip": False,
            "file_name": "hg38.cadd_snv.tsv.gz",
            "dir_name": "variants",
        },
    },
    GenomeVersion.CanFam3: {
        "reference_genome": {
            "url": "https://hgdownload.cse.ucsc.edu/goldenPath/canFam3/bigZips/canFam3.fa.gz",
            "file_type": FileType.FASTA,
            "gzip": True,
            "file_name": "canFam3.fasta",
            "dir_name": "genome",
        },
        "refgene_txt": {
            "url": "https://hgdownload.cse.ucsc.edu/goldenPath/canFam3/database/refGene.txt.gz",
            "file_type": FileType.TXT,
            "gzip": True,
            "file_name": "canfam3_refGene.txt",
            "dir_name": "genome",
        },
        "refgene_sql": {
            "url": "https://hgdownload.cse.ucsc.edu/goldenPath/canFam3/database/refGene.sql",
            "file_type": FileType.TXT,
            "gzip": False,
            "file_name": "canfam3_refGene.sql",
            "dir_name": "genome",
        },
        "genome_chrom_size": {
            "url": "https://hgdownload.cse.ucsc.edu/goldenPath/canFam3/bigZips/canFam3.chrom.sizes",
            "file_type": FileType.TXT,
            "gzip": False,
            "file_name": "canfam3.chrom.sizes",
            "dir_name": "genome",
        },
    },
}
