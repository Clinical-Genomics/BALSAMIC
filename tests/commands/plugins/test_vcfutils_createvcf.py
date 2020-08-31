from datetime import date

from BALSAMIC.commands.plugins.vcfutils import vcfheader
from BALSAMIC.commands.plugins.vcfutils import collect_vcf_info
from BALSAMIC.commands.plugins.vcfutils import collect_ref_info
from BALSAMIC.commands.plugins.vcfutils import readinput
from BALSAMIC.commands.plugins.vcfutils import createvcf


import re

def input_file():
    return "tests/test_data/vcf_tables/test_input.txt"

def vcf_output():
    return "tests/test_data/vcf_tables/test_createVCF_output.vcf.gz"

def reference_file():
    return "tests/test_data/vcf_tables/test_reference.vcf.gz"


def test_vcfheader_return_string():
    """test vcfheader for properly returning a VCF header"""
    #GIVEN current datetime
    current_time = date.today().strftime("%Y%m%d")
    valid_date_in_vcf = "##fileDate=" + current_time

    #WHEN calling vcfheader
    built_vcf_header = vcfheader()

    #THEN it should return a VCF header with a valid current date
    assert valid_date_in_vcf in built_vcf_header


def test_ensids_return_string():
    """ test ensembl ids in a reference vcf file """
    #GIVEN the ensembl ID
    info = 'GENE=AKT1_ENST00000555528'
    valid_ens_id = re.sub(r'(.*)(_ENST\d+)', r'\1', info)

    #WHEN substitute the ensembl_ids
    ens_id = collect_vcf_info(valid_ens_id)

    #THEN it should return a valid matching
    assert valid_ens_id in ens_id


def test_refinfo_return_list():
    """ test if the input file is a tab-delimited text file """
    #GIVEN the variant info fiels in input file
    valid_variant = "0.00119999998;SNP;p.Glu17Lys"
    allele_freq, variant_type, aa_hgvs = valid_variant.split(';')
    valid_info_variant = "VARIANT_TYPE=SNP;AA_HGVS=p.Glu17Lys;AF=1e-05"   
 
    #WHEN calling the ref file
    info = collect_ref_info(valid_variant)

    #THEN it should return a list of info fields
    assert valid_info_variant in info
 

