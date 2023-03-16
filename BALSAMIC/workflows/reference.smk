# syntax=python tabstop=4 expandtab
# coding: utf-8

import os
import logging
from pathlib import Path

from copy import deepcopy

from BALSAMIC.utils.rule import get_script_path
from BALSAMIC.utils.rule import get_reference_output_files
from BALSAMIC.utils.models import ReferenceMeta
from BALSAMIC.constants.reference import REFERENCE_FILES as REFERENCE_MODEL
from BALSAMIC.utils.cli import get_md5
from BALSAMIC.utils.cli import create_md5


LOG = logging.getLogger(__name__)

# explictly check if cluster_config dict has zero keys.
if len(cluster_config.keys()) == 0:
    cluster_config = config


# backward compatible genome version extraction from config
if 'genome_version' in config:
    genome_ver = config['genome_version']
else:
    genome_ver = 'hg19'

# essential path reference files
basedir = os.path.join(config['output'])
genome_dir = os.path.join(basedir, "genome")
vcf_dir = os.path.join(basedir, "variants")
vep_dir = os.path.join(basedir, "vep")
cosmicdb_key = config['cosmic_key']

# Set temporary dir environment variable
os.environ['TMPDIR'] = basedir

# indexable VCF files
# For future reference, if you delete this line pydantic fails in tests
# Don't know why, but don't delete it and keep deepcopy /A&H
REFERENCE_FILES = deepcopy(REFERENCE_MODEL)
indexable_vcf_files = get_reference_output_files(REFERENCE_FILES[genome_ver],
                                                 file_type='vcf',
                                                 gzip = True)

# intialize reference files
REFERENCE_FILES[genome_ver]['basedir'] = basedir
reference_file_model = ReferenceMeta.parse_obj(REFERENCE_FILES[genome_ver])

reference_genome_url = reference_file_model.reference_genome
dbsnp_url = reference_file_model.dbsnp
hc_vcf_1kg_url = reference_file_model.hc_vcf_1kg
mills_1kg_url = reference_file_model.mills_1kg
known_indel_1kg_url = reference_file_model.known_indel_1kg
vcf_1kg_url = reference_file_model.vcf_1kg
gnomad_url = reference_file_model.gnomad_variant
gnomad_tbi_url = reference_file_model.gnomad_variant_index
cosmicdb_url = reference_file_model.cosmicdb
wgs_calling_url = reference_file_model.wgs_calling
genome_chrom_size_url = reference_file_model.genome_chrom_size
refgene_txt_url = reference_file_model.refgene_txt
refgene_sql_url = reference_file_model.refgene_sql
rankscore_url = reference_file_model.rankscore
access_regions_url = reference_file_model.access_regions
delly_exclusion_url = reference_file_model.delly_exclusion
delly_mappability_url = reference_file_model.delly_mappability
delly_mappability_gindex_url = reference_file_model.delly_mappability_gindex
delly_mappability_findex_url = reference_file_model.delly_mappability_findex
ascat_gccorrection_url = reference_file_model.ascat_gccorrection
ascat_chryloci_url = reference_file_model.ascat_chryloci
clinvar_url = reference_file_model.clinvar
somalier_sites_url = reference_file_model.somalier_sites

# add secrets from config to items that need them
cosmicdb_url.secret=config['cosmic_key']

check_md5 = os.path.join(basedir, "reference.json.md5")

shell.executable("/bin/bash")
shell.prefix("set -eo pipefail; ")

singularity_image_path = config['singularity']['image_path']
singularity_images = [Path(singularity_image_path, image_name + ".sif").as_posix() for image_name in config["singularity"]["containers"].keys()]

##########################################################
# Generating Reference files for BALSAMIC pipeline
# Writing reference json file
##########################################################

rule all:
    input:
        singularity_images,
        reference_genome = reference_genome_url.get_output_file,
        bwa_index = expand(reference_genome_url.get_output_file + "{ext}", ext=['.amb','.ann','.bwt','.pac','.sa']),
        refgenome_fai = reference_genome_url.get_output_file + ".fai",
        refgenome_dict = reference_genome_url.get_output_file.replace("fasta","dict"),
        refseq_bed = refgene_txt_url.get_output_file.replace("txt", "flat") + ".bed",
        refseq_flat = refgene_txt_url.get_output_file.replace("txt", "flat"),
        refgene = refgene_txt_url.get_output_file,
        dbsnp_vcf = dbsnp_url.get_output_file + ".gz",
        th_genome_vcf = vcf_1kg_url.get_output_file + ".gz",
        tg_high_vcf = hc_vcf_1kg_url.get_output_file+ ".gz",
        mills_1kg = mills_1kg_url.get_output_file + ".gz",
        known_indel_1kg = known_indel_1kg_url.get_output_file + ".gz",
        gnomad_variant_vcf = gnomad_url.get_output_file,
        gnomad_variant_index = gnomad_tbi_url.get_output_file,
        cosmic_vcf = cosmicdb_url.get_output_file + ".gz",
        variants_idx = expand(os.path.join(vcf_dir,"{vcf}.gz.tbi"), vcf=indexable_vcf_files),
        vep = directory(vep_dir),
        wgs_calling = wgs_calling_url.get_output_file,
        genome_chrom_size = genome_chrom_size_url.get_output_file,
        rankscore = rankscore_url.get_output_file,
        access_regions = access_regions_url.get_output_file,
        delly_exclusion = delly_exclusion_url.get_output_file,
        delly_exclusion_converted = delly_exclusion_url.get_output_file.replace(".tsv", "_converted.tsv"),
        delly_mappability= delly_mappability_url.get_output_file,
        delly_mappability_gindex= delly_mappability_gindex_url.get_output_file,
        delly_mappability_findex= delly_mappability_findex_url.get_output_file,
        ascat_gccorrection = ascat_gccorrection_url.get_output_file,
        ascat_chryloci = ascat_chryloci_url.get_output_file,
        clinvar = clinvar_url.get_output_file + ".gz",
        somalier_sites = somalier_sites_url.get_output_file + ".gz",
    output:
        finished = os.path.join(basedir,"reference.finished"),
        reference_json = os.path.join(basedir, "reference.json"),
        check_md5 = check_md5
    log:
        os.path.join(basedir, "reference.json.log")
    run:
        import json
        from datetime import datetime

        today = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

        ref_json = dict()
        ref_json['reference'] = {
            "reference_genome": input.reference_genome,
            "dbsnp": input.dbsnp_vcf,
            "1kg_snps_all": input.th_genome_vcf,
            "1kg_snps_high": input.tg_high_vcf,
            "1kg_known_indel": input.known_indel_1kg,
            "mills_1kg": input.mills_1kg,
            "gnomad_variant": input.gnomad_variant_vcf,
            "cosmic": input.cosmic_vcf,
            "exon_bed": input.refseq_bed,
            "refflat": input.refseq_flat,
            "refGene": input.refgene,
            "wgs_calling_interval": input.wgs_calling,
            "genome_chrom_size": input.genome_chrom_size,
            "vep": input.vep,
            "rankscore": input.rankscore,
            "access_regions": input.access_regions,
            "delly_exclusion" : input.delly_exclusion,
            "delly_exclusion_converted" : input.delly_exclusion_converted,
            "delly_mappability": input.delly_mappability,
            "ascat_gccorrection" : input.ascat_gccorrection,
            "ascat_chryloci" : input.ascat_chryloci,
            "clinvar": input.clinvar,
            "somalier_sites": input.somalier_sites,
            "reference_access_date": today,
        }

        with open(str(output.reference_json), "w") as fh:
            json.dump(ref_json, fh, indent=4)

        create_md5(ref_json['reference'], output.check_md5)

        with open(str(output.finished), mode='w') as finish_file:
            finish_file.write('%s\n' % today )

###########################################################
# Download all singularity container images from dockerhub
###########################################################
wildcard_constraints:
    container_image = "|".join(list(config["singularity"]["containers"])),

def download_container_file(output_file: str):
    image_name = Path(output_file).stem
    docker_path = config["singularity"]["containers"][image_name]
    cmd = "singularity pull {}/{}.sif {}".format(config["singularity"]["image_path"],image_name,docker_path)
    shell(cmd)

rule download_container:
    output: Path(singularity_image_path, "{container_image}" + ".sif").as_posix(),
    run:
      download_container_file(output_file=output[0])

##########################################################
# Download the reference genome, variant db
##########################################################
download_content = [reference_genome_url, dbsnp_url, hc_vcf_1kg_url,
                    mills_1kg_url, known_indel_1kg_url, vcf_1kg_url,
                    wgs_calling_url, genome_chrom_size_url,
                    gnomad_url, gnomad_tbi_url,
                    cosmicdb_url, refgene_txt_url, refgene_sql_url, rankscore_url, access_regions_url,
                    delly_exclusion_url, delly_mappability_url, delly_mappability_gindex_url,
                    delly_mappability_findex_url, ascat_gccorrection_url, ascat_chryloci_url, clinvar_url,
                    somalier_sites_url]

download_dict = dict([(ref.get_output_file, ref) for ref in download_content])

def download_reference_file(output_file: str):
    import requests

    ref = download_dict[output_file]
    log_file = output_file + ".log"

    if ref.url.scheme == "gs":
        cmd = "export TMPDIR=/tmp; gsutil cp -L {} {} -".format(log_file,ref.url)
    else:
        cmd = "wget -a {} -O - {}".format(log_file,ref.url)

    if ref.secret:
        try:
            response = requests.get(ref.url,headers={'Authorization': 'Basic %s' % ref.secret})
            download_url = response.json()["url"]
        except:
            LOG.error("Unable to download {}".format(ref.url))
            raise
        cmd = "curl -o - '{}'".format(download_url)

    if ref.gzip:
        cmd += " | gunzip "

    cmd += " > {}".format(output_file)
    shell(cmd)
    ref.write_md5

ref_subdirs = set([ref.output_path for ref in download_content])
ref_files = set([ref.output_file for ref in download_content])

wildcard_constraints:
    ref_subdir="|".join(ref_subdirs),
    ref_file = "|".join(ref_files),


rule download_reference:
    output:
        Path("{ref_subdir}","{ref_file}").as_posix(),
    run:
        download_reference_file(output_file=output[0])




##########################################################
# Preprocess refseq file by fetching relevant columns and
# standardize the chr column
##########################################################

rule prepare_refgene:
    input:
        singularity_images,
        refgene_txt = refgene_txt_url.get_output_file,
        refgene_sql = refgene_sql_url.get_output_file,
        accessible_regions  = access_regions_url.get_output_file,
    params:
        refgene_sql_awk = get_script_path('refseq_sql.awk'),
    output:
        refflat = refgene_txt_url.get_output_file.replace("txt", "flat"),
        bed = refgene_txt_url.get_output_file.replace("txt", "flat") + ".bed",
    log:
        refgene_sql = os.path.join(basedir, "genome", "refgene_sql.log"),
        refgene_txt = os.path.join(basedir, "genome", "refgene_txt.log")
    singularity: Path(singularity_image_path, config["bioinfo_tools"].get("bedtools") + ".sif").as_posix()
    shell:
        """
header=$(awk -f {params.refgene_sql_awk} {input.refgene_sql});
(echo \"$header\"; cat {input.refgene_txt};) \
| csvcut -t -c chrom,exonStarts,exonEnds,name,score,strand,exonCount,txStart,txEnd,name2 \
| csvformat -T \
| bedtools expand -c 2,3 \
| awk '$1~/chr[1-9]/ && $1!~/[_]/' | cut -c 4- | sort -k1,1 -k2,2n > {output.bed};

awk -v OFS=\"\\t\" '$3!~/_/ {{ gsub(\"chr\",\"\",$3); $1=$13; print }}' {input.refgene_txt} \
| cut -f 1-11 > {output.refflat};
sed -i 's/chr//g' {input.refgene_txt};
sed -i 's/chr//g' {input.accessible_regions};
        """

##########################################################
# bgzip and tabix the vcf files that are vcf
##########################################################

rule bgzip_tabix:
    input:
        singularity_img = singularity_images,
        vcf = os.path.join(vcf_dir, "{vcf}.vcf")
    params:
        type = 'vcf',
    output:
        os.path.join(vcf_dir, "{vcf}.vcf.gz"),
        os.path.join(vcf_dir, "{vcf}.vcf.gz.tbi")
    log:
        os.path.join(vcf_dir, "{vcf}.vcf.gz_tbi.log")
    singularity: Path(singularity_image_path, config["bioinfo_tools"].get("tabix") + ".sif").as_posix()
    shell:
        """
bgzip {input.vcf} && tabix -p {params.type} {input.vcf}.gz 2> {log};
        """


##########################################################
# Create BWA Index for reference genome
##########################################################

rule bwa_index:
    input:
        singularity_img = singularity_images,
        reference_genome = reference_genome_url.get_output_file
    output:
        expand(reference_genome_url.get_output_file + "{ext}", ext=['.amb','.ann','.bwt','.pac','.sa'])
    log:
        reference_genome_url.get_output_file + ".bwa_index.log"
    singularity: Path(singularity_image_path, config["bioinfo_tools"].get("bwa") + ".sif").as_posix()
    shell:
        """
bwa index -a bwtsw {input.reference_genome} 2> {log};
        """

##########################################################
# Create index for fasta file - .fai
##########################################################

rule samtools_index_fasta:
    input:
        singularity_img = singularity_images,
        reference_genome = reference_genome_url.get_output_file
    output:
        reference_genome_url.get_output_file + ".fai"
    log:
        reference_genome_url.get_output_file + ".faidx.log"
    singularity: Path(singularity_image_path, config["bioinfo_tools"].get("samtools") + ".sif").as_posix()
    shell:
        """
samtools faidx {input.reference_genome} 2> {log};
        """


##########################################################
# create reference dictionary using picard
##########################################################

rule picard_ref_dict:
    input:
        singularity_img = singularity_images,
        reference_genome = reference_genome_url.get_output_file
    output:
        reference_genome_url.get_output_file.replace("fasta","dict")
    log:
        reference_genome_url.get_output_file + ".ref_dict.log"
    singularity: Path(singularity_image_path, config["bioinfo_tools"].get("picard") + ".sif").as_posix()
    shell:
        """
picard CreateSequenceDictionary REFERENCE={input.reference_genome} OUTPUT={output} 2> {log};
        """


##########################################################
# ENSEMBL VEP - download and install vep package,
#                 cache conversion
##########################################################

rule vep_install:
    input:
        singularity_img = singularity_images
    params:
        species = "homo_sapiens_merged",
        assembly = "GRCh37" if genome_ver == 'hg19' else "GRCh38",
        plugins = "all",
    output:
        directory(vep_dir)
    log:
        os.path.join(vep_dir, "vep_install_cache.log")
    singularity: Path(singularity_image_path, config["bioinfo_tools"].get("ensembl-vep") + ".sif").as_posix()
    shell:
        """
vep_install --SPECIES {params.species} \
--AUTO cfp \
--ASSEMBLY {params.assembly} \
--CACHEDIR {output} \
--PLUGINS {params.plugins} \
--NO_HTSLIB --CONVERT --NO_UPDATE 2> {log}; 
        """

##########################################################
# Remove chr from delly exclusion
##########################################################

rule prepare_delly_exclusion:
    input:
        singularity_img = singularity_images,
        delly_exclusion = delly_exclusion_url.get_output_file,
    output:
        delly_exclusion_converted = delly_exclusion_url.get_output_file.replace(".tsv", "_converted.tsv"),
    log:
        os.path.join(basedir, "genome", "delly_exclusion.log"),
    singularity: Path(singularity_image_path, config["bioinfo_tools"].get("delly") + ".sif").as_posix()
    shell:
        """
sed 's/chr//g' {input.delly_exclusion} > {output.delly_exclusion_converted} 2> {log}
        """
