# syntax=python tabstop=4 expandtab
# coding: utf-8

import os
import hashlib
import logging

from datetime import date

from BALSAMIC.utils.rule import get_script_path
from BALSAMIC.utils.rule import get_reference_output_files 
from BALSAMIC.utils.models import ReferenceMeta
from BALSAMIC.utils.constants import REFERENCE_FILES 

LOG = logging.getLogger(__name__)

# explictly check if cluster_config dict has zero keys.
if len(cluster_config.keys()) == 0:
    cluster_config = config

current_day = date.today()

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
ascat_gccorrection_url = reference_file_model.ascat_gccorrection
ascat_chryloci_url = reference_file_model.ascat_chryloci

# add secrets from config to items that need them
cosmicdb_url.secret=config['cosmic_key']

check_md5 = os.path.join(basedir, "reference.json.md5")

shell.executable("/bin/bash")
shell.prefix("set -eo pipefail; ")

def get_md5(filename):
    hash_md5 = hashlib.md5()
    with open(str(filename), 'rb') as fh:
        for chunk in iter(lambda: fh.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def create_md5(reference, check_md5):
    """ create a md5 file for all reference data"""
    with open(check_md5, 'w') as fh:
        for key, value in reference.items():
            if os.path.isfile(value):
                fh.write( get_md5(value) + ' ' + value + '\n')

singularity_image = config['singularity']['image']

##########################################################
# Generating Reference files for BALSAMIC pipeline
# Writing reference json file 
##########################################################

rule all:
    input:
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
        ascat_gccorrection = ascat_gccorrection_url.get_output_file,
        ascat_chryloci = ascat_chryloci_url.get_output_file
    output:
        finished = os.path.join(basedir,"reference.finished"),
        reference_json = os.path.join(basedir, "reference.json"),
        check_md5 = check_md5
    log:
        os.path.join(basedir, "reference.json.log")
    run:
        import json

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
            "ascat_gccorrection" : input.ascat_gccorrection,
            "ascat_chryloci" : input.ascat_chryloci
        }

        with open(str(output.reference_json), "w") as fh:
            json.dump(ref_json, fh, indent=4)
        
        create_md5(ref_json['reference'], output.check_md5)

        shell("date +'%Y-%M-%d T%T %:z' > {output.finished}") 


##########################################################
# Download the reference genome, variant db 
##########################################################
download_content = [reference_genome_url, dbsnp_url, hc_vcf_1kg_url,
                    mills_1kg_url, known_indel_1kg_url, vcf_1kg_url,
                    wgs_calling_url, genome_chrom_size_url,
                    gnomad_url, gnomad_tbi_url,
                    cosmicdb_url, refgene_txt_url, refgene_sql_url, rankscore_url, access_regions_url,
                    delly_exclusion_url, ascat_gccorrection_url, ascat_chryloci_url]

rule download_reference:
    output:
        expand("{output}", output=[ref.get_output_file for ref in download_content])
    run:
        import requests

        for ref in download_content:
            output_file = ref.get_output_file
            log_file = output_file + ".log"

            if ref.url.scheme == "gs":
                cmd = "export TMPDIR=/tmp; gsutil cp -L {} {} -".format(log_file, ref.url)
            else:
                cmd = "wget -a {} -O - {}".format(log_file, ref.url)

            if ref.secret:
                try:
                    response = requests.get(ref.url, headers={'Authorization': 'Basic %s' % ref.secret })
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


rule prepare_refgene:
    input:
        refgene_txt = refgene_txt_url.get_output_file,
        refgene_sql = refgene_sql_url.get_output_file,
        accessible_regions  = access_regions_url.get_output_file,
    params:
        refgene_sql_awk = get_script_path('refseq_sql.awk'),
        conda_env = config["bioinfo_tools"].get("bedtools")
    output:
        refflat = refgene_txt_url.get_output_file.replace("txt", "flat"),
        bed = refgene_txt_url.get_output_file.replace("txt", "flat") + ".bed",
    log:
        refgene_sql = os.path.join(basedir, "genome", "refgene_sql.log"),
        refgene_txt = os.path.join(basedir, "genome", "refgene_txt.log")
    singularity: Path(singularity_image, config["bioinfo_tools"].get("bedtools") + ".sif").as_posix() 
    shell:
        """
source activate {params.conda_env};
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
        os.path.join(vcf_dir, "{vcf}.vcf")
    params:
        type = 'vcf',
        conda_env = config["bioinfo_tools"].get("tabix")
    output:
        os.path.join(vcf_dir, "{vcf}.vcf.gz"),
        os.path.join(vcf_dir, "{vcf}.vcf.gz.tbi")
    log:
        os.path.join(vcf_dir, "{vcf}.vcf.gz_tbi.log")
    singularity: Path(singularity_image, config["bioinfo_tools"].get("tabix") + ".sif").as_posix() 
    shell:
        """
source activate {params.conda_env};
bgzip {input} && tabix -p {params.type} {input}.gz 2> {log};
        """


##########################################################
# Create BWA Index for reference genome
##########################################################

rule bwa_index:
    input:
        reference_genome_url.get_output_file
    params:
        conda_env = config["bioinfo_tools"].get("bwa")
    output:
        expand(reference_genome_url.get_output_file + "{ext}", ext=['.amb','.ann','.bwt','.pac','.sa'])
    log:
        reference_genome_url.get_output_file + ".bwa_index.log"
    singularity: Path(singularity_image, config["bioinfo_tools"].get("bwa") + ".sif").as_posix() 
    shell:
        """
source activate {params.conda_env};
bwa index -a bwtsw {input} 2> {log};
        """

##########################################################
# Create index for fasta file - .fai
##########################################################

rule samtools_index_fasta:
    input:
        reference_genome_url.get_output_file
    params:
        conda_env = config["bioinfo_tools"].get("samtools")
    output:
        reference_genome_url.get_output_file + ".fai"
    log:
        reference_genome_url.get_output_file + ".faidx.log"
    singularity: Path(singularity_image, config["bioinfo_tools"].get("samtools") + ".sif").as_posix() 
    shell:
        """
source activate {params.conda_env};
samtools faidx {input} 2> {log};
        """


##########################################################
# create reference dictionary using picard
# 
##########################################################

rule picard_ref_dict:
    input:
        reference_genome_url.get_output_file
    params:
        conda_env = config["bioinfo_tools"].get("picard")
    output:
        reference_genome_url.get_output_file.replace("fasta","dict")
    log:
        reference_genome_url.get_output_file + ".ref_dict.log"
    singularity: Path(singularity_image, config["bioinfo_tools"].get("picard") + ".sif").as_posix() 
    shell:
        """
source activate {params.conda_env};
picard CreateSequenceDictionary REFERENCE={input} OUTPUT={output} 2> {log};
        """


##########################################################
# ENSEMBL VEP - download and install vep package, 
#                 cache coversion
##########################################################

rule vep_install:
    params:
        species = "homo_sapiens_merged",
        assembly = "GRCh37" if genome_ver == 'hg19' else "GRCh38",
        plugins = "all",
        conda_env = config["bioinfo_tools"].get("ensembl-vep")
    output:
        directory(vep_dir)
    log:
        os.path.join(vep_dir, "vep_install_cache.log")
    singularity: Path(singularity_image, config["bioinfo_tools"].get("ensembl-vep") + ".sif").as_posix() 
    shell:
        """
source activate {params.conda_env};
vep_install --SPECIES {params.species} \
--AUTO cfp \
--ASSEMBLY {params.assembly} \
--CACHEDIR {output} \
--PLUGINS {params.plugins} \
--NO_HTSLIB --CONVERT --NO_UPDATE 2> {log}; 
        """


rule prepare_delly_exclusion:
    input:
        delly_exclusion = delly_exclusion_url.get_output_file,
    output:
        delly_exclusion_converted = delly_exclusion_url.get_output_file.replace(".tsv", "_converted.tsv"),
    log:
        os.path.join(basedir, "genome", "delly_exclusion.log"),
    singularity:
        Path(singularity_image, config["bioinfo_tools"].get("delly") + ".sif").as_posix()
    shell:
        """
sed 's/chr//g' {input.delly_exclusion} > {output.delly_exclusion_converted} 2> {log}
        """
