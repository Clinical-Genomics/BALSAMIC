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

genome_ver = config['genome_version']

# essential path reference files
basedir = os.path.join(config['output'])
genome_dir = os.path.join(basedir, "genome")

# Set temporary dir environment variable
os.environ['TMPDIR'] = basedir

REFERENCE_FILES = deepcopy(REFERENCE_MODEL)

# intialize reference files
REFERENCE_FILES[genome_ver]['basedir'] = basedir
reference_file_model = ReferenceMeta.parse_obj(REFERENCE_FILES[genome_ver])
reference_genome_url = reference_file_model.reference_genome
genome_chrom_size_url = reference_file_model.genome_chrom_size
refgene_txt_url = reference_file_model.refgene_txt
refgene_sql_url = reference_file_model.refgene_sql

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
        genome_chrom_size = genome_chrom_size_url.get_output_file,
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
            "exon_bed": input.refseq_bed,
            "refflat": input.refseq_flat,
            "refGene": input.refgene,
            "genome_chrom_size": input.genome_chrom_size,
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

rule download_container:
    output: singularity_images
    run:
      for image_name, docker_path in config["singularity"]["containers"].items():
          cmd = "singularity pull {}/{}.sif {}".format(config["singularity"]["image_path"], image_name, docker_path)
	  shell(cmd)

##########################################################
# Download the reference genome, variant db
##########################################################
download_content = [reference_genome_url,  genome_chrom_size_url, refgene_txt_url, refgene_sql_url]

rule download_reference:
    output:
        expand("{output}", output=[ref.get_output_file for ref in download_content])
    run:
        import requests

        for ref in download_content:
            output_file = ref.get_output_file
            log_file = output_file + ".log"

            cmd = "wget -a {} -O - {}".format(log_file, ref.url)

            if ref.gzip:
                cmd += " | gunzip "

            cmd += " > {}".format(output_file)
            shell(cmd)
            ref.write_md5

##########################################################
# Preprocess refseq file by fetching relevant columns and
# standardize the chr column
##########################################################

rule prepare_refgene:
    input:
        singularity_images,
        refgene_txt = refgene_txt_url.get_output_file,
        refgene_sql = refgene_sql_url.get_output_file,
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
| awk '$1~/chr[1-9]/ && $1!~/[_]/' | sort -k1,1 -k2,2n > {output.bed};

awk -v OFS=\"\\t\" '$3!~/_/ {{ gsub(\"chr\",\"chr\",$3); $1=$13; print }}' {input.refgene_txt} \
| cut -f 1-11 > {output.refflat};
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

