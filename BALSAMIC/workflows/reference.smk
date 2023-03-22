# syntax=python tabstop=4 expandtab
# coding: utf-8

import logging
import os
from pathlib import Path

from BALSAMIC.constants.cache import FileType
from BALSAMIC.constants.paths import BALSAMIC_DIR
from BALSAMIC.constants.workflow_rules import SNAKEMAKE_RULES
from BALSAMIC.models.cache_models import CacheConfigModel
from BALSAMIC.utils.io import write_finish_file
from BALSAMIC.utils.rule import get_threads

LOG = logging.getLogger(__name__)

# Balsamic cache configuration model
cache_config: CacheConfigModel = CacheConfigModel.parse_obj(config)

# Temporary directory environment
os.environ["TMPDIR"] = cache_config.references_dir.as_posix()
shell.executable("/bin/bash")
shell.prefix("set -eo pipefail; ")

# Rules to include
for rule in SNAKEMAKE_RULES["cache"]:
    include: Path(BALSAMIC_DIR, rule).as_posix()
LOG.info(f"The rules {SNAKEMAKE_RULES['cache']} will be included in the reference workflow")

rule all:
    """Target rule for Balsamic cache generation."""
    input:
        expand(
            Path(cache_config.containers_dir, "{singularity_image}." + FileType.SIF).as_posix(),
            singularity_image=cache_config.containers.keys(),
        ),
        #cache_config.get_reference_paths(),
    output:
        finish_file=Path(cache_config.references_dir, "reference.finish"),
    threads: get_threads(cluster_config, "all")
    run:
        write_finish_file(output.finish_file)


# ##########################################################
# # Download the reference genome, variant db
# ##########################################################
# download_content = [
#     reference_genome_url,
#     dbsnp_url,
#     hc_vcf_1kg_url,
#     mills_1kg_url,
#     known_indel_1kg_url,
#     vcf_1kg_url,
#     wgs_calling_url,
#     genome_chrom_size_url,
#     gnomad_url,
#     gnomad_tbi_url,
#     cosmicdb_url,
#     refgene_txt_url,
#     refgene_sql_url,
#     rankscore_url,
#     access_regions_url,
#     delly_exclusion_url,
#     delly_mappability_url,
#     delly_mappability_gindex_url,
#     delly_mappability_findex_url,
#     ascat_gccorrection_url,
#     ascat_chryloci_url,
#     clinvar_url,
#     somalier_sites_url,
# ]
#
# download_dict = dict([(ref.get_output_file, ref) for ref in download_content])
#
#
# def download_reference_file(output_file: str):
#     import requests
#
#     ref = download_dict[output_file]
#     log_file = output_file + ".log"
#
#     if ref.url.scheme == "gs":
#         cmd = "export TMPDIR=/tmp; gsutil cp -L {} {} -".format(log_file, ref.url)
#     else:
#         cmd = "wget -a {} -O - {}".format(log_file, ref.url)
#
#     if ref.secret:
#         try:
#             response = requests.get(
#                 ref.url, headers={"Authorization": "Basic %s" % ref.secret}
#             )
#             download_url = response.json()["url"]
#         except:
#             LOG.error("Unable to download {}".format(ref.url))
#             raise
#         cmd = "curl -o - '{}'".format(download_url)
#
#     if ref.gzip:
#         cmd += " | gunzip "
#
#     cmd += " > {}".format(output_file)
#     shell(cmd)
#     ref.write_md5
#
#
# ref_subdirs = set([ref.output_path for ref in download_content])
# ref_files = set([ref.output_file for ref in download_content])
#
#
# wildcard_constraints:
#     ref_subdir="|".join(ref_subdirs),
#     ref_file="|".join(ref_files),
#
#
# rule download_reference:
#     output:
#         Path("{ref_subdir}", "{ref_file}").as_posix(),
#     run:
#         download_reference_file(output_file=output[0])
#
#
# ##########################################################
# # Preprocess refseq file by fetching relevant columns and
# # standardize the chr column
# ##########################################################
#
#
# rule prepare_refgene:
#     input:
#         cache_config.container_paths,
#         refgene_txt=refgene_txt_url.get_output_file,
#         refgene_sql=refgene_sql_url.get_output_file,
#         accessible_regions=access_regions_url.get_output_file,
#     params:
#         refgene_sql_awk=get_script_path("refseq_sql.awk"),
#     output:
#         refflat=refgene_txt_url.get_output_file.replace("txt", "flat"),
#         bed=refgene_txt_url.get_output_file.replace("txt", "flat") + ".bed",
#     log:
#         refgene_sql=Path(reference_dir, "genome", "refgene_sql.log"),
#         refgene_txt=Path(reference_dir, "genome", "refgene_txt.log"),
#     singularity:
#         Path(
#             singularity_image_path, config["bioinfo_tools"].get("bedtools") + ".sif"
#         ).as_posix()
#     shell:
#         """
# header=$(awk -f {params.refgene_sql_awk} {input.refgene_sql});
# (echo \"$header\"; cat {input.refgene_txt};) \
# | csvcut -t -c chrom,exonStarts,exonEnds,name,score,strand,exonCount,txStart,txEnd,name2 \
# | csvformat -T \
# | bedtools expand -c 2,3 \
# | awk '$1~/chr[1-9]/ && $1!~/[_]/' | cut -c 4- | sort -k1,1 -k2,2n > {output.bed};
#
# awk -v OFS=\"\\t\" '$3!~/_/ {{ gsub(\"chr\",\"\",$3); $1=$13; print }}' {input.refgene_txt} \
# | cut -f 1-11 > {output.refflat};
# sed -i 's/chr//g' {input.refgene_txt};
# sed -i 's/chr//g' {input.accessible_regions};
#         """
#
#
# ##########################################################
# # bgzip and tabix the vcf files that are vcf
# ##########################################################
#
#
# rule bgzip_tabix:
#     input:
#         singularity_img=cache_config.container_paths,
#         vcf=Path(vcf_dir, "{vcf}.vcf"),
#     params:
#         type="vcf",
#     output:
#         Path(vcf_dir, "{vcf}.vcf.gz"),
#         Path(vcf_dir, "{vcf}.vcf.gz.tbi"),
#     log:
#         Path(vcf_dir, "{vcf}.vcf.gz_tbi.log"),
#     singularity:
#         Path(
#             singularity_image_path, config["bioinfo_tools"].get("tabix") + ".sif"
#         ).as_posix()
#     shell:
#         """
# bgzip {input.vcf} && tabix -p {params.type} {input.vcf}.gz 2> {log};
#         """
#
#
# ##########################################################
# # Create BWA Index for reference genome
# ##########################################################
#
#
# rule bwa_index:
#     input:
#         singularity_img=cache_config.container_paths,
#         reference_genome=reference_genome_url.get_output_file,
#     output:
#         expand(
#             reference_genome_url.get_output_file + "{ext}",
#             ext=[".amb", ".ann", ".bwt", ".pac", ".sa"],
#         ),
#     log:
#         reference_genome_url.get_output_file + ".bwa_index.log",
#     singularity:
#         Path(
#             singularity_image_path, config["bioinfo_tools"].get("bwa") + ".sif"
#         ).as_posix()
#     shell:
#         """
# bwa index -a bwtsw {input.reference_genome} 2> {log};
#         """
#
#
# ##########################################################
# # Create index for fasta file - .fai
# ##########################################################
#
#
# rule samtools_index_fasta:
#     input:
#         singularity_img=cache_config.container_paths,
#         reference_genome=reference_genome_url.get_output_file,
#     output:
#         reference_genome_url.get_output_file + ".fai",
#     log:
#         reference_genome_url.get_output_file + ".faidx.log",
#     singularity:
#         Path(
#             singularity_image_path, config["bioinfo_tools"].get("samtools") + ".sif"
#         ).as_posix()
#     shell:
#         """
# samtools faidx {input.reference_genome} 2> {log};
#         """
#
#
# ##########################################################
# # create reference dictionary using picard
# ##########################################################
#
#
# rule picard_ref_dict:
#     input:
#         singularity_img=cache_config.container_paths,
#         reference_genome=reference_genome_url.get_output_file,
#     output:
#         reference_genome_url.get_output_file.replace("fasta", "dict"),
#     log:
#         reference_genome_url.get_output_file + ".ref_dict.log",
#     singularity:
#         Path(
#             singularity_image_path, config["bioinfo_tools"].get("picard") + ".sif"
#         ).as_posix()
#     shell:
#         """
# picard CreateSequenceDictionary REFERENCE={input.reference_genome} OUTPUT={output} 2> {log};
#         """
#
#
# ##########################################################
# # ENSEMBL VEP - download and install vep package,
# #                 cache conversion
# ##########################################################
#
#
# rule vep_install:
#     input:
#         singularity_img=cache_config.container_paths,
#     params:
#         species="homo_sapiens_merged",
#         assembly="GRCh37" if genome_version == "hg19" else "GRCh38",
#         plugins="all",
#     output:
#         directory(vep_dir),
#     log:
#         Path(vep_dir, "vep_install_cache.log"),
#     singularity:
#         Path(
#             singularity_image_path, config["bioinfo_tools"].get("ensembl-vep") + ".sif"
#         ).as_posix()
#     shell:
#         """
# vep_install --SPECIES {params.species} \
# --AUTO cfp \
# --ASSEMBLY {params.assembly} \
# --CACHEDIR {output} \
# --PLUGINS {params.plugins} \
# --NO_HTSLIB --CONVERT --NO_UPDATE 2> {log};
#         """
#
#
# ##########################################################
# # Remove chr from delly exclusion
# ##########################################################
#
#
# rule prepare_delly_exclusion:
#     input:
#         singularity_img=cache_config.container_paths,
#         delly_exclusion=delly_exclusion_url.get_output_file,
#     output:
#         delly_exclusion_converted=delly_exclusion_url.get_output_file.replace(
#             ".tsv", "_converted.tsv"
#         ),
#     log:
#         Path(reference_dir, "genome", "delly_exclusion.log"),
#     singularity:
#         Path(
#             singularity_image_path, config["bioinfo_tools"].get("delly") + ".sif"
#         ).as_posix()
#     shell:
#         """
# sed 's/chr//g' {input.delly_exclusion} > {output.delly_exclusion_converted} 2> {log}
#         """
