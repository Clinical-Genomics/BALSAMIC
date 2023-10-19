
from BALSAMIC.assets.scripts.preprocess_gens import (
    extract_variant_info,
    get_valid_variants,
    extract_coverage_line_values
)

import logging
from typing import Dict
from io import StringIO

def test_extract_variant_info():
    """test extraction of variant information from a vcf line from DNAscope."""

    # GIVEN VALID VARIANT
    valid_variant = "1\t100\trs1\tT\tC\t389.77\t.\tINFO\tGT:AD:DP:GQ:PL\t0/1:9,14:23:99:418,0,257"
    variant_info = extract_variant_info(valid_variant)

    # THEN AF should be found and correctly calculated
    assert variant_info["af"] > 0.6 and variant_info["af"] < 0.61

    # GIVEN INVALID VARIANT
    invalid_variant = "1\t200\t.\tCAAA\tCAAAA,C\t0.00\tLowQual\tINFO\tGT:AD:DP:GQ:PL\t0/0:0,0,0:0:0:0,0,0,3,3,19"
    variant_info = extract_variant_info(invalid_variant)
    assert variant_info == None

def test_get_valid_variants(caplog):
    """test extraction of variant information from a vcf line from DNAscope."""

    # GIVEN LIST OF VARIANT STRINGS
    variant_list = ["1\t100\trs1\tT\tC\t389.77\t.\tINFO\tGT:AD:DP:GQ:PL\t0/1:9,14:23:99:418,0,257",
                    "1\t200\t.\tCAAA\tCAAAA,C\t0.00\tLowQual\tINFO\tGT:AD:DP:GQ:PL\t0/0:0,0,0:0:0:0,0,0,3,3,19",
                    "25\t100\trs1\tT\tC\t389.77\t.\tINFO\tGT:AD:DP:GQ:PL\t0/1:9,14:23:99:418,0,257"]

    # Call the function
    variant_dict: Dict = get_valid_variants(variant_list)

    # Check if we get the expected variant dict
    assert variant_dict[0] == {"chr": "1", "start": "100", "ref": "T", "alt": "C",
                               "sample": "0/1:9,14:23:99:418,0,257", "af": 0.608696}

    # Check if any warning messages were logged
    assert any(record.levelname == "WARNING" for record in caplog.records)
    # Check if a specific warning message was logged
    assert any("Can't calc AF for a number of variants: 1." in record.message for record in caplog.records)
    assert any("A number of variants have illegal chromosomes and will be skipped: {'25': 1}." in record.message for record in caplog.records)


def test_extract_coverage_line_values():
    """test extraction coverage region values from a coverage line.."""

    # GIVEN A COVERAGE LINE
    coverage_line = "1\t13401\t13500\t-0.956941"
    chrom, start, end, log2_ratio = extract_coverage_line_values(coverage_line)
    assert chrom == "1"
    assert start == 13401
    assert end == 13500
    assert log2_ratio == -0.956941