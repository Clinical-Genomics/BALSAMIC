from datetime import date

from BALSAMIC.commands.plugins.vcfutils import vcfheader
from BALSAMIC.commands.plugins.vcfutils import collect_vcf_info
from BALSAMIC.commands.plugins.vcfutils import collect_ref_info
from BALSAMIC.commands.plugins.vcfutils import readinput
from BALSAMIC.commands.plugins.vcfutils import createvcf
from click.testing import CliRunner

import re
import pytest


def test_readinput_return_list(input_file):
    """ test input file for properly returning the required fields """
    #GIVEN input file path
    valid_input_info = {
        'COSV62571334:AKT1:p.E17K': '0.118;SNP;p.Glu17Lys',
        'COSV51765161:EGFR:p.L858R': '0.106;SNP;p.Leu858Arg',
        'COSV51765492:EGFR:p.T790M': '0.116;SNP;p.Thr790Met',
        'COSV56056643:BRAF:p.V600E': '0.104;SNP;p.Val600Glu'
    }

    #WHEN calling readinput
    build_read_input = readinput(input_file)

    #THEN it should return a input info with dict value
    assert valid_input_info == build_read_input


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
    """ test fields in reference file """
    #GIVEN the variant info fiels in input file
    valid_variant = "0.00119999998;SNP;p.Glu17Lys"
    allele_freq, variant_type, aa_hgvs = valid_variant.split(';')
    valid_info_variant = "VARIANT_TYPE=SNP;AA_HGVS=p.Glu17Lys;AF=1e-05"

    #WHEN calling the ref file
    info = collect_ref_info(valid_variant)

    #THEN it should return a list of info fields
    assert valid_info_variant in info


def test_createvcf(input_file, reference_file, output_file):
    """ test createvcf """
    runner = CliRunner()
    result = runner.invoke(
        createvcf, ["-i", input_file, "-r", reference_file, "-o", output_file])
    assert result.exit_code == 0, "all files exits"
