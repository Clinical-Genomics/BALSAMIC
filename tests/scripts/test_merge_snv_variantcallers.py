import pytest
from BALSAMIC.assets.scripts.merge_snv_variantcallers import (
    parse_chromosome,
    merge_filters,
    merge_info_fields,
    merge_headers,
    VCFHeaderMergeError,
    merge_variants,
)
from typing import Tuple, List
import os


def test_merge_variants(test_vardict_vcf: str, test_tnscope_vcf: str):
    """Tests correct merging of variants from VCF files which share 1 variant and have 1 unique each."""

    # GIVEN test VarDict and TNscope VCF containing 2 variants each, sharing 1

    # WHEN merging variants that have the same chromosome, pos, ref and alt
    merged_variants: List[str] = merge_variants(test_vardict_vcf, test_tnscope_vcf)

    # THEN a total of 3 variants should exist in the final list of variants
    assert len(merged_variants) == 3


def test_merge_variantheaders(test_vardict_vcf: str, test_tnscope_vcf: str):
    """Tests correct merging headers from VarDict and TNscope VCFs"""

    # GIVEN test VarDict and TNscope VCF

    # WHEN merging headers
    merged_header: List[str] = merge_headers(test_vardict_vcf, test_tnscope_vcf)

    # THEN rows with different descriptions should be merged from both VCFs
    assert (
        '##FORMAT=<ID=AF,Number=.,Type=Float,Description="vcf1: Allele Frequency | vcf2: Allele fraction of the event in the tumor">'
        in merged_header
    )

    # THEN only 1 chromosome entry should exist per chromosome
    chrom1_rows = [
        hdr_row
        for hdr_row in merged_header
        if "contig=<ID=1,length=249250621" in hdr_row
    ]
    assert len(chrom1_rows) == 1

    # THEN unique VarDict headers should exist in the merged VCF
    assert (
        '##INFO=<ID=VD,Number=1,Type=Integer,Description="Variant Depth">'
        in merged_header
    )

    # THEN unique TNscope headers should exist in the merged VCF
    assert (
        '##INFO=<ID=NLOD,Number=1,Type=Float,Description="Normal LOD score">'
        in merged_header
    )

    # THEN new header row should exist in the merged VCF
    vcf1_name = os.path.basename(test_vardict_vcf)
    vcf2_name = os.path.basename(test_tnscope_vcf)
    assert (
        f"##merge_snv_variantcallers=values in merged INFO fields are listed in the order of the input files: first from {vcf1_name}, then from {vcf2_name}"
        in merged_header
    )


def test_merge_headers_mismatch(
    test_tnscope_vcf: str, test_vardict_vcf_non_matching_header: str
):
    """Tests correct error capture when merging headers with different numbers of samples."""
    with pytest.raises(
        VCFHeaderMergeError, match="Error: Variant headers in .* do not match"
    ):
        merge_headers(test_tnscope_vcf, test_vardict_vcf_non_matching_header)


def test_sorting_and_parsing_of_chromosome():
    """Tests accurate sorting order of non-numbered chromosome."""

    variants1 = {
        ("X", "100", "A", "T"): {},  # unique to v 1
        ("2", "100", "A", "T"): {},  # should be merged
        ("1", "100", "A", "T"): {},  # unique to v 1
        ("Y", "100", "A", "T"): {},  # should be merged
        ("1", "400", "A", "T"): {},  # should be merged
        ("M", "100", "A", "T"): {},  # should be merged
        ("3", "100", "A", "T"): {},  # unique to v 1
    }
    variants2 = {
        ("X", "105", "A", "T"): {},  # unique to v 2
        ("2", "100", "A", "T"): {},  # should be merged
        ("1", "150", "A", "T"): {},  # unique to v 2
        ("Y", "100", "A", "T"): {},  # should be merged
        ("1", "400", "A", "T"): {},  # should be merged
        ("3", "450", "A", "T"): {},  # unique to v 2
        ("M", "100", "A", "T"): {},  # should be merged
        ("3", "101", "A", "C"): {},  # unique to v 2
    }

    all_keys: List[Tuple[str, str, str, str]] = sorted(
        set(variants1.keys()).union(variants2.keys()),
        key=lambda x: (parse_chromosome(x[0]), int(x[1])),
    )

    expected_sorted_keys = [
        ("1", "100", "A", "T"),
        ("1", "150", "A", "T"),
        ("1", "400", "A", "T"),
        ("2", "100", "A", "T"),
        ("3", "100", "A", "T"),
        ("3", "101", "A", "C"),
        ("3", "450", "A", "T"),
        ("X", "100", "A", "T"),
        ("X", "105", "A", "T"),
        ("Y", "100", "A", "T"),
        ("M", "100", "A", "T"),
    ]

    assert all_keys == expected_sorted_keys


def test_merge_conflicting_filters():
    """Tests that filters are correctly merged when conflicting."""

    # GIVEN conflicting filters for the same variant
    filter1 = "PASS"
    filter2 = "in_normal;germline_risk"

    # WHEN merging filters
    merged_filters = merge_filters(filter1, filter2)

    # THEN PASS filter should be removed
    assert merged_filters == "germline_risk;in_normal"


def test_merge_filters():
    """Tests that filters are correctly merged when agreeing."""

    # GIVEN conflicting filters for the same variant
    filter1 = "PASS"
    filter2 = "PASS"

    # WHEN merging filters
    merged_filters = merge_filters(filter1, filter2)

    # THEN only one PASS should remain
    assert merged_filters == "PASS"


def test_merge_info_fields():
    """Test that all types of info field mergings are done as expected."""

    # Test merging unique key-value pairs
    assert merge_info_fields(["DP=10", "AF=0.5"]) == "DP=10;AF=0.5"

    # Test merging duplicate keys with different values
    assert merge_info_fields(["DP=10", "DP=20"]) == "DP=10;DP_LIST=10,20"

    # Test handling of key-only fields
    assert merge_info_fields(["P0.01Likely", "DP=10"]) == "P0.01Likely;DP=10"

    # Test merging key-only fields with duplicate keys
    assert merge_info_fields(["PASS", "PASS", "DP=10"]) == "PASS;DP=10"

    # Test merging multiple values with a mix of key-only and key-value pairs
    assert (
        merge_info_fields(["NLODF=126.72", "DP=10", "DP=15", "PASS"])
        == "NLODF=126.72;DP=10;PASS;DP_LIST=10,15"
    )

    # Test empty input
    assert merge_info_fields([]) == ""

    # Test handling of mixed formatting
    assert (
        merge_info_fields(["AF=0.1", "AF=0.2", "FILTER", "FILTER"])
        == "AF=0.1;FILTER;AF_LIST=0.1,0.2"
    )

    # Test handling of multiple unique key-value pairs
    assert merge_info_fields(["MQ=60", "DP=10", "AF=0.3"]) == "MQ=60;DP=10;AF=0.3"

    # Test case with multiple duplicates
    assert (
        merge_info_fields(
            ["DP=10", "DP=20", "DP=30", "FOUND_IN=vardict", "FOUND_IN=tnscope"]
        )
        == "DP=10;FOUND_IN=vardict,tnscope;DP_LIST=10,20,30"
    )
