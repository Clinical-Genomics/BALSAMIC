# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import click
import re
from datetime import date
from datetime import datetime
from cyvcf2 import VCF

@click.group()
def VCFutils():
	""" Commands to process VCF files """
	
@VCFutils.command()
@click.option('-i', '--input_file', required=True, type=click.STRING,help='tab-seperated reference text file')
@click.option('--print/--no-print', default=False, help='tab-seperated reference text file')
@click.pass_context

 
def readinput(ctx, input_file, print):
    
    """ Input text file processing """
    """ 
    Tab-delimited input file should contain the following columns: 
    ['Gene_ID', 'COSMIC_ID', 'Mutation_ID', 'Variant_type', 'AA. Change', 'AA_HGVS', 'Average_AF%']
    Information is stored as dictionary where; Keys: Mutation_ID:Gene_ID:AA_change Values: Average_AF%;AA_HGVS;Variant_type 
    """    

    AF = {}
    ip_file = open(input_file,'r')
    header = ip_file.readline().strip().split('\t')
    for ln in iter(ip_file):
        ln = ln.strip().split('\t')
        ref_keys = ln[header.index('Mutation_ID')] + ':' + ln[header.index('Gene_ID')] + ':' + ln[header.index('AA_Change')]
        AF[ref_keys] = ln[header.index('Average_AF%')] + ';' + ln[header.index('Variant_type')] + ';'+ ln[header.index('AA_HGVS')]
    ip_file.close()
    if print:
        click.echo(AF)
    else:
        return AF



def vcfheader(output_file):

    """ Write standard header lines for VCF output """
    flDate = date.today().strftime("%Y%m%d")
    vcf_op = open(output_file,'w')
    vcf_op.write("##fileformat=VCFv4.2\n"
    "##fileDate="+flDate+"\n"
    "##source=NA\n"
    "##reference=NA\n"
    "##contig=<ID=20,length=NA,assembly=NA,md5=NA,species=\"Homo sapiens\",taxonomy=x>\n"
    "##INFO=<ID=GENE,Number=1,Type=String,Description=\"Gene_name\">\n"
    "##INFO=<ID=AA,Number=1,Type=String,Description=\"Aminoacid mutation\">\n"
    "##INFO=<ID=VARIANT_TYPE,Number=1,Type=String,Description=\"Type of mutation\">\n"
    "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Average Allele frequency fraction\">\n"
    "##INFO=<ID=AA_HGVS,Number=1,Type=Flag,Description=\"AA_HGVS\">\n"
    "##INFO=<ID=COSMIC,Number=1,Type=String,Description=\"Genomic Mutation Identifier\">\n"
    "##INFO=<ID=LEGACY_ID,Number=1,Type=String,Description=\"Legacy Mutation Identifier\">\n"
    "##INFO=<ID=CDS,Number=1,Type=String,Description=\"CDS mutation\">\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    return(vcf_op)



@VCFutils.command()
@click.option('-i', '--input_file', required=True, type=click.STRING,help='tab-seperated reference text file')
@click.option('-r', '--reference_file', required=True, type=click.STRING,help='cosmic database file')
@click.option('-o', '--output_file', required=True, type=click.STRING,help='Output file name')
@click.pass_context

def createvcf(ctx, input_file,reference_file, output_file):

    """ Filter input variants from reference VCF """
    """ Exact information which is collected as keys in read_input() is retrieved and matched"""
    """ Compare the matching keys in both files and extract additional information """
    start_time = datetime.now()
    vcf_op = vcfheader(output_file)
    #AF = read_input(input_file)
    AF = ctx.invoke(readinput,input_file=input_file)
    filtered_variants = []
    for variant in VCF(reference_file):
        """ skipping transctipt_ids """
        gene_symbol = re.sub(r'(.*)(_ENST.*)',r'\1', variant.INFO.get('GENE'))
        vcf_id = variant.ID + ':' + gene_symbol + ':' + variant.INFO.get('AA')
        vcf_id_value = AF.get(vcf_id)
        if vcf_id_value:
            #click.echo('key is:' + str(vcf_id) +';' +'value is:'+ str(vcf_id_value))
            af, VARIANT_TYPE, AA_HGVS = vcf_id_value.split(';')
            variant_info = str(variant).split('\t')
            """ Ensembl transcript IDs and CNTs values are not required; remove it"""
            """ Convert af values to fraction and join all relevant information """
            variant_info = [re.sub(r'(.*)(_ENST\d+;)',r'\1;', i) for i in variant_info]
            variant_info = [re.sub(r'(.*)(;CNT=\d+\n)',r'\1', i) for i in variant_info]
            add_info = ";".join(["VARIANT_TYPE="+ str(VARIANT_TYPE), "AA_HGVS="+str(AA_HGVS),"AF="+str(round(float(af)/100,5))])
            info = "\t".join(variant_info) + ';' + add_info
            if info not in filtered_variants:
               	#click.echo(info)
                filtered_variants.append(info)
                vcf_op.write(info + '\n')
    click.echo("vcf file created")
    end_time = datetime.now()
    click.echo("Runtime:" + '{}'.format(end_time - start_time))
    return filtered_variants

@VCFutils.command()
@click.option('-f1', '--input_file1', required=True, type=click.STRING,help='Reference VCF file')
@click.option('-f2', '--input_file2', required=True, type=click.STRING,help='VEP annotated VCF file')
@click.option('-o', '--output_file', required=True, type=click.STRING,help='Output file name')
@click.pass_context

def compareVCF(ctx,input_file1,input_file2,output_file):
    """ Compare two VCF files (truthset vs validate set) """
