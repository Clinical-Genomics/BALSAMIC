# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
""" Command-line utilities for processing VCF files """

import re
from datetime import date
from datetime import datetime
import click
from cyvcf2 import VCF


@click.group()
def vcfutils():
    """ Commands to process VCF files """


##Tab-delimited input file containing specific header format.Information is stored as dictionary.
def readinput(text_file):
    """ Input text file processing. Outputs a dictionary """
    input_info = {}
    with open(text_file, 'r') as input_file:
        header = input_file.readline().strip().split('\t')
        for lines in iter(input_file):
            lines = lines.strip().split('\t')
            ref_keys = ":".join([
                lines[header.index('Mutation_ID')],
                lines[header.index('Gene_ID')],
                lines[header.index('AA_Change')]
            ])
            input_info[ref_keys] = ";".join([
                lines[header.index('Average_AF%')],
                lines[header.index('Variant_type')],
                lines[header.index('AA_HGVS')]
            ])
    return input_info


def vcfheader():
    """ Write standard header lines for VCF output """
    file_date = date.today().strftime("%Y%m%d")
    vcf_header = "{}".format("""##fileformat=VCFv4.2
##fileDate=""" + file_date + """
##source=NA
##reference=NA
##contig=<ID=20,length=NA,assembly=NA,md5=NA,species=\"Homo sapiens\",taxonomy=x>
##INFO=<ID=GENE,Number=1,Type=String,Description=\"Gene_name\">
##INFO=<ID=AA,Number=1,Type=String,Description=\"Aminoacid mutation\">
##INFO=<ID=VARIANT_TYPE,Number=1,Type=String,Description=\"Type of mutation\">
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Average Allele frequency fraction\">
##INFO=<ID=AA_HGVS,Number=1,Type=Flag,Description=\"AA_HGVS\">
##INFO=<ID=COSMIC,Number=1,Type=String,Description=\"Genomic Mutation Identifier\">
##INFO=<ID=LEGACY_ID,Number=1,Type=String,Description=\"Legacy Mutation Identifier\">
##INFO=<ID=CDS,Number=1,Type=String,Description=\"CDS mutation\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO""")
    return vcf_header


def collect_ref_info(variant):
    """ Collect necessary variant information from input file"""
    allele_freq, variant_type, aa_hgvs = variant.split(';')
    info = ";".join([
        "VARIANT_TYPE=" + str(variant_type), "AA_HGVS=" + str(aa_hgvs),
        "AF=" + str(round(float(allele_freq) / 100, 5))
    ])
    return info


def collect_vcf_info(variant):
    """ Clean reference vcf and collect info from required fields """
    info = str(variant).split('\t')
    info = [re.sub(r'(.*)(_ENST\d+;)', r'\1;', i) for i in info]
    info = [re.sub(r'(.*)(;CNT=\d+\n)', r'\1', i) for i in info]
    return info


@vcfutils.command()
@click.option('-i',
              '--input_file',
              required=True,
              type=click.Path(exists=True),
              help='tab-seperated reference text file')
@click.option('-r',
              '--reference_file',
              required=True,
              type=click.Path(exists=True),
              help='cosmic database file')
@click.option('-o',
              '--output_file',
              required=True,
              type=click.Path(exists=True),
              help='Output file name')
def createvcf(input_file, reference_file, output_file):
    """ Filter input variants from reference VCF """
    ## Exact information which is collected as keys in read_input() is retrieved
    ## Compare the matching keys in both files and extract additional information

    start_time = datetime.now()
    allele_freq = readinput(input_file)
    filtered_variants = []
    with open(output_file, 'w') as vcf_output:
        vcf_output.write(vcfheader())
        for variant in VCF(reference_file):
            gene_symbol = re.sub(r'(.*)(_ENST.*)', r'\1',
                                 variant.INFO.get('GENE'))
            vcf_id = ":".join(
                [variant.ID, gene_symbol,
                 variant.INFO.get('AA')])
            vcf_id_value = allele_freq.get(vcf_id)
            if vcf_id_value:
                variant_info = collect_vcf_info(str(variant))
                reference_info = collect_ref_info(vcf_id_value)
                info = "\t".join(variant_info) + ';' + reference_info
                if info not in filtered_variants:
                    filtered_variants.append(info)
                    vcf_output.write(info + '\n')
    end_time = datetime.now()
    click.echo("VCF file created. Total Runtime:" +
               '{}'.format(end_time - start_time))
    return filtered_variants


@vcfutils.command()
@click.option('-f1',
              '--input_file1',
              required=True,
              type=click.Path(exists=True),
              help='Reference VCF file')
@click.option('-f2',
              '--input_file2',
              required=True,
              type=click.Path(exists=True),
              help='VEP annotated VCF file')
@click.option('-o',
              '--output_file',
              required=True,
              type=click.Path(exists=True),
              help='Output file name')
def comparevcf(input_file1, input_file2, output_file):
    """ Compare two VCF files (truthset vs validate set) """
    if input_file1:
        pass
    if input_file2:
        pass
    if output_file:
        pass
