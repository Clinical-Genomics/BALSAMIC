from datetime import date

from BALSAMIC.commands.plugins.vcfutils import vcfheader

def test_vcfheader_return_string():
    """test vcfheader for properly returning a VCF header"""
    #GIVEN current datatime
    current_time = date.today().strftime("%Y%m%d")
    valid_date_in_vcf = "##fileDate=" + current_time

    #WHEN calling vcfheader
    built_vcf_header = vcfheader()

    #THEN it should return a VCF header with a valid current date
    assert valid_date_in_vcf in built_vcf_header
