from BALSAMIC.assets.scripts.preprocess_gens import (
    extract_variant_info,
    get_valid_variants,
    extract_coverage_line_values,
)

from typing import Dict


def test_extract_variant_info(valid_dnascope_variant, invalid_dnascope_variant_no_ad):
    """test extraction of variant information from a vcf line from DNAscope."""

    # GIVEN VALID VARIANT

    # WHEN extracting the variant information and calculating allele frequencies
    variant_info = extract_variant_info(valid_dnascope_variant)

    # THEN AF should be correctly calculated
    assert variant_info["af"] > 0.6 and variant_info["af"] < 0.61

    # GIVEN INVALID VARIANT

    # WHEN extracting the variant information and calculating allele frequencies
    variant_info = extract_variant_info(invalid_dnascope_variant_no_ad)

    # THEN no information should be returned
    assert variant_info == None


def test_get_valid_variants(
    caplog,
    valid_dnascope_variant,
    invalid_dnascope_variant_no_ad,
    invalid_dnascope_variant_illegal_chrom,
):
    """test extraction of variant information from a vcf line from DNAscope."""

    # GIVEN LIST OF 2 INVALID AND 1 VALID VARIANT STRINGS
    variant_list = [
        valid_dnascope_variant,
        invalid_dnascope_variant_no_ad,
        invalid_dnascope_variant_illegal_chrom,
    ]

    # WHEN extracting information for valid variants
    variant_dict: Dict = get_valid_variants(variant_list)

    # THEN the first variant should be correctly extracted into this structure
    assert variant_dict[0] == {
        "chr": "1",
        "start": "100",
        "ref": "T",
        "alt": "C",
        "sample": "0/1:9,14:23:99:418,0,257",
        "af": 0.608696,
    }
    # THEN the remaining variants should throw WARNINGS
    assert any(record.levelname == "WARNING" for record in caplog.records)
    # THEN one variant should fail calculation of allele frequency and be ignored
    assert any(
        "Can't calc AF for a number of variants: 1." in record.message
        for record in caplog.records
    )
    # THEN one variant should have an illegal chromosome name and be ignored
    assert any(
        "A number of variants have illegal chromosomes and will be skipped: {'25': 1}."
        in record.message
        for record in caplog.records
    )


def test_extract_coverage_line_values():
    """test extraction coverage region values from a coverage line.."""

    # GIVEN A COVERAGE LINE
    coverage_line = "1\t13401\t13500\t-0.956941"

    # WHEN extracting information from string
    chrom, start, end, log2_ratio = extract_coverage_line_values(coverage_line)

    # THEN the coverage region variables should have the expected values
    assert chrom == "1"
    assert start == 13401
    assert end == 13500
    assert log2_ratio == -0.956941
