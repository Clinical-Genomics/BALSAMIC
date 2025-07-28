from BALSAMIC.assets.scripts.whitelist_snvs import (
    parse_info,
    partial_match,
    onc_clnsig,
    matches_whitelist,
    whitelist_variants,
)
import pytest
from io import StringIO


# --- parse_info ---
def test_parse_info():
    info = "DP=100;AF=0.5;CLNSIG=Pathogenic;OTHER=Yes"
    expected = {"DP": "100", "AF": "0.5", "CLNSIG": "Pathogenic", "OTHER": "Yes"}
    assert parse_info(info) == expected


def test_parse_info_empty():
    assert parse_info("") == {}


# --- partial_match ---
@pytest.mark.parametrize(
    "query,target,expected",
    [
        ("A", "A", True),
        ("A", "B", False),
        ("...T", "GCT", True),
        ("C...", "CAT", True),
        ("...T...", "AGCTG", True),
    ],
)
def test_partial_match(query, target, expected):
    assert partial_match(query, target) == expected


# --- onc_clnsig ---
def test_onc_clnsig_positive():
    info = "CLNSIG=Pathogenic;DP=50"
    assert onc_clnsig(info) is True

    info = "ONC=Oncogenic;AF=0.1"
    assert onc_clnsig(info) is True


def test_onc_clnsig_negative():
    info = "CLNSIG=Benign;DP=10"
    assert onc_clnsig(info) is False

    info = "DP=100;AF=0.5"
    assert onc_clnsig(info) is False


# --- matches_whitelist ---
def test_matches_whitelist_exact_match():
    variant = {"chrom": "1", "pos": "100", "ref": "A", "alt": "T"}
    white = {"chrom": "1", "pos": "100", "ref": "A", "alt": "T"}
    assert matches_whitelist(variant, white)


def test_matches_whitelist_partial_match():
    variant = {"chrom": "1", "pos": "200", "ref": "GCT", "alt": "A"}
    white = {"chrom": "1", "pos": "200", "ref": "...T", "alt": "A"}
    assert matches_whitelist(variant, white)


def test_matches_whitelist_mismatch():
    variant = {"chrom": "2", "pos": "100", "ref": "G", "alt": "C"}
    white = {"chrom": "1", "pos": "100", "ref": "G", "alt": "C"}
    assert not matches_whitelist(variant, white)


# --- whitelist_variants ---
def test_whitelist_variants_onc_clnsig():
    variants = {
        ("1", "100", "A", "T"): {
            "chrom": "1",
            "pos": "100",
            "id": ".",
            "ref": "A",
            "alt": "T",
            "qual": ".",
            "filter": "LowQual",
            "info": "BaseQRankSum=-0.000;ClippingRankSum=-0.000;CLNSIG=Pathogenic;DP=10;ExcessHet=3.0103",
            "format": "GT",
            "samples": ["0/1"],
        }
    }
    result = whitelist_variants(variants.copy())
    v = result[("1", "100", "A", "T")]
    assert "WhitelistStatus=ClinvarPathogenicOncogenic" in v["info"]
    assert "WhitelistedFilters=LowQual" in v["info"]
    assert v["filter"] == "PASS"


def test_whitelist_variants_with_list_match():
    variants = {
        ("1", "123", "G", "A"): {
            "chrom": "1",
            "pos": "123",
            "id": ".",
            "ref": "G",
            "alt": "A",
            "qual": ".",
            "filter": "TriallelicSite",
            "info": "DP=30",
            "format": "GT",
            "samples": ["0/1"],
        }
    }
    whitelist = {1: {"chrom": "1", "pos": "123", "ref": "G", "alt": "A"}}
    result = whitelist_variants(variants.copy(), whitelist)
    v = result[("1", "123", "G", "A")]
    assert "WhitelistStatus=ClinicalList" in v["info"]
    assert "WhitelistedFilters=TriallelicSite" in v["info"]
    assert v["filter"] == "PASS"
