# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

shell.prefix('set -eo pipefail;')

# *************************************************************
# Workflow to perform UMI tag extraction using 
# Sentieon tools and calling somatic variants 
# Analysis run in 3 modes -
#   a) UMI structure defined on both reads (ds)
#   b) UMI structure defined on first read (ss_r1)
#   c) UMI structure defined on second read (ss_r2)
# *************************************************************
# Workflow references: 
# https://support.sentieon.com/appnotes/umi/  
# *************************************************************


configfile: "config.yaml"
localrules: all

import os
import glob
from os.path import join 

base_dir = config['base_dir']
wdir = join(base_dir, config['case_ID'], 'UMI_run')
workdir: wdir


try:
    SENTIEON_LICENSE = os.environ["SENTIEON_LICENSE"]
    SENTIEON_INSTALL_DIR = os.environ["SENTIEON_INSTALL_DIR"]
except Expection as error:
    print('ERROR: Set environment variable to run Sentieon')
    raise
    
fastq_dir = join(base_dir, config['case_ID'], 'fastq' )
(SAMPLES_PE, READS,) = glob_wildcards(os.path.join(fastq_dir,"{sample}_R_{read}.fastq.gz"))
ip_reads = [join(fastq_dir, "{sample}_R_1.fastq.gz"), join(fastq_dir, "{sample}_R_2.fastq.gz")]
VAR_CALLER  = ['TNscope', 'vardict', 'mutect2']
MODE = ['ds', 'ss_r1', 'ss_r2']

# define outputs for rules
singularity_image = config['singularity_img']

tnscope_vcf = expand('TNscope/{sample}_TNscope_{mode}.vcf.gz', sample=SAMPLES_PE, mode=MODE)

AF_tables = expand('{var_caller}/tables/{sample}_AFtable_{var_caller}_{mode}.txt', sample=SAMPLES_PE, mode=MODE, var_caller = VAR_CALLER)

vardict_vcf = expand('vardict/{sample}_vardict_{mode}.vcf.gz', sample=SAMPLES_PE, mode=MODE)

mutect2_vcf = expand('mutect2/{sample}_mutect2_{mode}.vcf.gz mutect2/{sample}_mutect2_contamination_table_{mode}.vcf.gz'.split(),sample=SAMPLES_PE, mode=MODE)

vep_annotate = expand('vep/{sample}_{var_caller}_vep_{mode}.vcf.gz', sample=SAMPLES_PE, var_caller=VAR_CALLER, mode=MODE)


rule all:
    input: tnscope_vcf + AF_tables + vardict_vcf + vep_annotate + mutect2_vcf
 
# Extract UMI-tags
rule umi_extract:
    input:
        reads = ip_reads
    output:
        ds_umi = 'umi/{sample}_umiextract_ds.fastq.gz',
        ss_r1_umi = 'umi/{sample}_umiextract_ss_r1.fastq.gz',
        ss_r2_umi = 'umi/{sample}_umiextract_ss_r2.fastq.gz'
    params:
        sentieon_exec = SENTIEON_INSTALL_DIR + 'bin/sentieon',
        sentieon_lic = SENTIEON_LICENSE,
        ds_params= ['-d', '3M2S+T,3M2S+T'],
        ss_r1_params= [' ', '3M2S+T,5S+T'],
        ss_r2_params = [' ', '5S+T,3M2S+T'],
        sample_id = '{sample}'
    log: 'logs/{sample}_umiextract.log'
    benchmark:  'benchmarks/{sample}_umiextract.benchmark'
    message: "UMI tag extraction using Sentieon for sample {params.sample_id}"
    shell:
        "export SENTIEON_LICENSE={params.sentieon_lic}\n"
        "{params.sentieon_exec} umi extract {params.ds_params} {input.reads} -o {output.ds_umi} &> {log} \n"
        "{params.sentieon_exec} umi extract {params.ss_r1_params} {input.reads} -o {output.ss_r1_umi} &>> {log} \n"
        "{params.sentieon_exec} umi extract {params.ss_r2_params} {input.reads} -o {output.ss_r2_umi} &>> {log} \n"

# Align the UMI-extracted reads
rule umi_align:
    input:
        ref_fa = config['references']['genomefa'],
        fastq_umi = 'umi/{sample}_umiextract_{mode}.fastq.gz'
    output:
        align_umi = 'umi/{sample}_umialign_{mode}.sam'
    params:
        sentieon_exec = SENTIEON_INSTALL_DIR + 'bin/sentieon',
        sentieon_lic = SENTIEON_LICENSE,
        sheader = "'@RG\\tID:Group\\tSM:{sample}\\tLB:TargetPanel\\tPL:ILLUMINA'",
        ip_bases = '1000000',
        sample_id = '{sample}',
        mode = '{mode}'
    threads: 16
    log: 'logs/{sample}_umialign_{mode}.log'
    benchmark : 'benchmarks/{sample}_umialign_{mode}.benchmark'
    message: "Aligning of UMI extracted reads with bwa mem, sorting for sample {params.sample_id} {params.mode}"
    shell:
        "export SENTIEON_LICENSE={params.sentieon_lic}\n"
        "{params.sentieon_exec} bwa mem -R {params.sheader} -K {params.ip_bases} -p -t {threads} -C {input.ref_fa} {input.fastq_umi} > {output.align_umi} 2> {log}\n"


# UMI-consensus calling
rule umi_consensus:
    input:
        sam_consensus = 'umi/{sample}_umialign_{mode}.sam',
    output:
        fastq_consensus = 'consensus/{sample}_consensuscall_{mode}.fastq.gz',
    params:
        sentieon_exec = SENTIEON_INSTALL_DIR + 'bin/sentieon',
        sentieon_lic = SENTIEON_LICENSE,
        tag = 'XR',
        ip_format = 'SAM',
        sample_id = '{sample}',
        mode = '{mode}'
    threads: 16
    log: 'logs/{sample}_consensuscall_{mode}.log'
    benchmark : 'benchmarks/{sample}_consensuscall_{mode}.benchmark'
    message: "Consensus molecule creation for sample {params.sample_id} {params.mode}"
    shell:
        "export SENTIEON_LICENSE={params.sentieon_lic}\n"
        "{params.sentieon_exec} umi consensus -t {threads} -i {input.sam_consensus} -o {output.fastq_consensus} --input_format {params.ip_format} --umi_tag {params.tag} --read_name_prefix 'UMI-' &> {log}\n"


# Alignment of consensus reads 
rule consensus_bam:
    input:
        ref_fa = config['references']['genomefa'],
        fq_consensus = 'consensus/{sample}_consensuscall_{mode}.fastq.gz',
    output:
        align_consensus = 'consensus/{sample}_consensusalign_{mode}.bam',
    params:
        sentieon_exec = SENTIEON_INSTALL_DIR + 'bin/sentieon',
        sentieon_lic = SENTIEON_LICENSE,
        sheader = "'@RG\\tID:Group\\tSM:{sample}\\tLB:TargetPanel\\tPL:ILLUMINA'",
        ip_bases = '1000000',
        sample_id = '{sample}',
        mode = '{mode}'
    threads: 16
    log: 'logs/{sample}_consensusalign_{mode}.log'
    benchmark : 'benchmarks/{sample}_consensusalign_{mode}.benchmark'
    message: "Mapping of consensus reads with the bwa mem, sorting for sample {params.sample_id} {params.mode}"
    shell:
        "export SENTIEON_LICENSE={params.sentieon_lic}\n"
        "{params.sentieon_exec} bwa mem -R {params.sheader} -t {threads} -K {params.ip_bases} -p -C {input.ref_fa} {input.fq_consensus} | {params.sentieon_exec} util sort -r {input.ref_fa} --sam2bam -o {output.align_consensus} -i - 2> {log}\n"


# Variant-calling using TNscope
rule tnscope:
    input:
        bam = 'consensus/{sample}_consensusalign_{mode}.bam',
        ref_fa = config['references']['genomefa'],
        bed = config['references']['targetPanel'],
        dbsnp = config['references']['dbsnp']
    output:
         vcf = 'TNscope/{sample}_TNscope_{mode}.vcf.gz',
    params:
        sentieon_exec = SENTIEON_INSTALL_DIR + 'bin/sentieon',
        sentieon_lic = SENTIEON_LICENSE,
        algo = 'TNscope',
        tAF = '0.0005',
        TL = '0.5',
        detect = 'sv',
        error_rate = '5',
        prune_factor = '3',
        sample_id = '{sample}',
        mode = '{mode}'
    threads: 16
    log: 'logs/{sample}_TNscope_{mode}.log'
    benchmark : 'benchmarks/{sample}_TNscope_{mode}.benchmark'
    message: "Calling SNVs using TNscope for sample {params.sample_id} {params.mode}"
    shell:
        "export SENTIEON_LICENSE={params.sentieon_lic}\n"
        "{params.sentieon_exec} driver -t {threads} -r {input.ref_fa} -i {input.bam} --algo {params.algo} --tumor_sample {params.sample_id} --dbsnp {input.dbsnp} --min_tumor_allele_frac {params.tAF} --filter_t_alt_frac {params.tAF} --min_init_tumor_lod {params.TL} --disable_detector {params.detect} --max_error_per_read {params.error_rate} --pcr_indel_model NONE --prune_factor {params.prune_factor} {output.vcf}  &> {log}\n"



# Variant-calling using vardict
rule vardict:
    input:
        bam = 'consensus/{sample}_consensusalign_{mode}.bam',
        ref_fa = config['references']['genomefa'],
        bed = config['references']['targetPanel']
    output:
         vardict = 'vardict/{sample}_vardict_{mode}.vcf.gz',
    params:
        conda = "BALSAMIC_py36" ,
        af = "0.0005",
        sample_id = '{sample}',
        mode = '{mode}',
        vardict = "-c 1 -S 2 -E 3 -g 4 -r 1 -F 0",
        var2vcf = '-E'    
    singularity: singularity_image
    threads: 8
    log: 
        'logs/{sample}_vardict_{mode}.log'
    benchmark:
        'benchmarks/{sample}_vardict_{mode}.benchmark'
    message:
        'Variant calling using Vardict for sample {params.sample_id} {params.mode}' 
    shell:
        "source activate {params.conda}\n"
        "vardict -G {input.ref_fa} -f {params.af} -N {params.sample_id} -b {input.bam} {params.vardict} {input.bed} | teststrandbias.R | var2vcf_valid.pl {params.var2vcf} -f {params.af} -N {params.sample_id} | bgzip > {output.vardict}\n"
        "tabix -p vcf {output.vardict}\n"
        "source deactivate"

# Variant-calling using mutect2
rule mutect2:
    input: 
        bam = 'consensus/{sample}_consensusalign_{mode}.bam',
        ref_fa = config['references']['genomefa'],
        bed = config['references']['targetPanel'],
        dbsnp = config['references']['dbsnp']
    output:
        temp_vcf = temp('mutect2/{sample}_mutect2_unfiltered_{mode}.vcf.gz'), 
        mutect2_vcf = 'mutect2/{sample}_mutect2_{mode}.vcf.gz',
        cont_tbl = 'mutect2/{sample}_mutect2_contamination_table_{mode}.vcf.gz',
    params:
        conda = 'D_UMI_APJ',
        sample_id = '{sample}',
        mode = '{mode}'
    threads: 8
    log: 
        'logs/{sample}_mutect2_{mode}.log'
    benchmark:
        'benchmarks/{sample}_mutect2_{mode}.benchmark'
    message:
        'Variant calling using mutect2 for sample {params.sample_id} {params.mode}'
    shell:
        "source activate {params.conda}\n"
        "gatk Mutect2 -R {input.ref_fa} -L {input.bed} -I {input.bam} --germline-resource {input.dbsnp} -O {output.temp_vcf}\n"
        "gatk FilterMutectCalls -V {output.temp_vcf} --contamination-table {output.cont_tbl} -O {output.mutect2_vcf}\n"  

# VEP annotations
rule vep_somatic:
    input:
        vcf = '{var_caller}/{sample}_{var_caller}_{mode}.vcf.gz',
    output:
        vep = 'vep/{sample}_{var_caller}_vep_{mode}.vcf.gz',
    params:
        conda = 'BALSAMIC_py36', 
        vep_cache = config['vep_cache'],
        default_options = '--compress_output bgzip --vcf --everything --allow_non_variant --dont_skip --buffer_size 10000 --format vcf --offline --variant_class --merged --cache --verbose --force_overwrite',
        sample_id = '{sample}',
        var_caller = '{var_caller}',
        mode = '{mode}'
    threads: 8 
    singularity: singularity_image
    log:
        'logs/{sample}_{var_caller}_{mode}_vep.log'
    benchmark:
        'benchmarks/{sample}_{var_caller}_{mode}_vep.benchmark'
    message:
        'Annotating {params.var_caller} VCF file with VEP for sample {params.sample_id} {params.mode}'
    shell:
        "source activate {params.conda}\n"
        "vep_path=$(dirname $(readlink -e $(which vep)))\n"
        "export PERL5LIB=\n"
        "vep --dir $vep_path --dir_cache {params.vep_cache} --dir_plugins $vep_path --input_file {input.vcf} --output_file {output.vep} --fork {threads} {params.default_options}\n"
        "tabix -p vcf -f {output.vep}\n"
        "source deactivate"

# Generate tables for AF scatterplots
rule calculate_AF:
    input:
        vcf = '{var_caller}/{sample}_{var_caller}_{mode}.vcf.gz',
    output: 
        AF = '{var_caller}/tables/{sample}_AFtable_{var_caller}_{mode}.txt',
    params:
        validated_set= config['background_snvs'],
        sample_id = '{sample}',
        mode='{mode}'
    threads: 4
    log:'logs/{sample}_AFcalculate_{var_caller}_{mode}.log'
    benchmark: 'benchmarks/{sample}_AFcalculate_{var_caller}_{mode}.benchmark'
    message: "Creating Allelic frequency table from VCF file for sample {params.sample_id} {params.mode}"
    shell:
        "bash /home/proj/development/cancer/UMI/utils/bcf2AF.sh {params.validated_set} {input.vcf} {params.sample_id} {output.AF} &> {log}\n"


