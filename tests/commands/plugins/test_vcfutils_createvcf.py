from datetime import date

from BALSAMIC.commands.plugins.vcfutils import vcfheader
from BALSAMIC.commands.plugins.vcfutils import collect_vcf_info
import re

def test_vcfheader_return_string():
    """test vcfheader for properly returning a VCF header"""
    #GIVEN current datatime
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

    #THEN it should return a valid matchi
    assert valid_ens_id in ens_id
