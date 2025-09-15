"""Tests for Balsamic variant filters classes."""
from math import isclose

import pytest
from typing import List

from BALSAMIC.constants.analysis import SequencingType, BioinfoTools, AnalysisType
from BALSAMIC.models.params import VCFFilter
from BALSAMIC.constants.variant_filters import (
    BaseSNVFilters,
    WgsSNVFilters,
    TgaSNVFilters,
    TgaUmiSNVFilters,
    get_tag_and_filtername,
)


def test_tga_quality_vardict_filters():
    """Test filter_criteria for excluding TNscope specific filters."""

    # GIVEN TGA SNV tumor only and excluding TNscope specific filters
    filters = TgaSNVFilters.filter_criteria(
        category="quality", variant_caller=BioinfoTools.VARDICT
    )

    # WHEN extracting filter names
    retrieved_filter_names = {f.filter_name for f in filters}

    # THEN vardict specific filter is in filters
    assert "balsamic_low_mq" in retrieved_filter_names
    # THEN TNscope specific filter is NOT in filters
    assert "balsamic_high_strand_oddsratio" not in retrieved_filter_names
    # THEN general quality filter is in filters
    assert "balsamic_low_tumor_dp" in retrieved_filter_names


def test_tga_quality_paired_filters():
    """Test filter_criteria for assigning Paired specific filters in TGA class."""

    # GIVEN TGA SNV paired filters
    filters = TgaSNVFilters.filter_criteria(
        category="quality", analysis_type=AnalysisType.PAIRED
    )

    # WHEN extracting filter names
    retrieved_filter_names = {f.filter_name for f in filters}

    # THEN matched normal filter is retrieved
    assert "in_normal" in retrieved_filter_names


def test_tga_quality_single_filters():
    """Test filter_criteria for assigning Single specific filters in TGA class."""

    # GIVEN TGA SNV tumor only filters
    filters = TgaSNVFilters.filter_criteria(
        category="quality", analysis_type=AnalysisType.SINGLE
    )

    # WHEN extracting filter names
    retrieved_filter_names = {f.filter_name for f in filters}

    # THEN matched normal filter is not retrieved
    assert "in_normal" not in retrieved_filter_names
    # THEN single specific filter is retrieved
    assert "balsamic_high_strand_oddsratio" in retrieved_filter_names


def test_wgs_quality_paired_filters():
    """Test filter_criteria for assigning Paired specific filters in WGS class."""

    # GIVEN WGS tumor normal matched filters
    filters = WgsSNVFilters.filter_criteria(
        category="quality", analysis_type=AnalysisType.PAIRED
    )

    # WHEN extracting filter names
    retrieved_filter_names = {f.filter_name for f in filters}

    # THEN matched normal filter is retrieved
    assert "in_normal" in retrieved_filter_names


def test_wgs_quality_single_filters():
    """Test filter_criteria for assigning Single specific filters in WGS class."""
    # GIVEN WGS tumor only filters
    filters = WgsSNVFilters.filter_criteria(
        category="quality", analysis_type=AnalysisType.SINGLE
    )

    # WHEN extracting filter names
    retrieved_filter_names = {f.filter_name for f in filters}

    # THEN matched normal filter is not retrieved
    assert "in_normal" not in retrieved_filter_names
    # THEN single specific filter is retrieved
    assert "balsamic_low_strand_read_counts" in retrieved_filter_names


def test_tga_umi_snv_filters():
    """Test correct filter retrieval in TGA UMI SNV filters class."""

    # GIVEN TGA UMI tumor normal matched filters
    filters = TgaUmiSNVFilters.filter_criteria(
        category="quality", analysis_type=AnalysisType.PAIRED
    )
    # WHEN extracting filter names
    retrieved_filter_names = {f.filter_name for f in filters}

    # THEN matched normal filter is retrieved
    assert "in_normal" in retrieved_filter_names

    # GIVEN TGA UMI research filters
    filters = TgaUmiSNVFilters.filter_criteria(category="research")

    # WHEN extracting filter names
    retrieved_filter_names = [f.filter_name for f in filters]

    # THEN all research filters are retrieved
    assert ["SWEGENAF", "balsamic_high_pop_freq"] == retrieved_filter_names


def test_wes_tga_read_depth_snv_filter():
    """Test correct retrieval of WES specific DP value"""

    # GIVEN WES and TGA filters
    wes_filters = TgaSNVFilters.get_filters(category="quality", exome=True)
    tga_filters = TgaSNVFilters.get_filters(category="quality", exome=False)

    # WHEN extracting DP specific filter
    wes_dp_filter = get_tag_and_filtername(wes_filters, "balsamic_low_tumor_dp")
    tga_dp_filter = get_tag_and_filtername(tga_filters, "balsamic_low_tumor_dp")

    # THEN DP filter name should be correct
    assert wes_dp_filter[1] == "balsamic_low_tumor_dp"
    # THEN DP threshold for WES should be 20
    assert isclose(wes_dp_filter[0], 20)
    # THEN DP threshold for TGA should be 50
    assert isclose(tga_dp_filter[0], 50)


def test_get_bcftools_filter_string_soft_normals_exclusion():
    """Test soft_filter_normals excludes matched normal filters from filter string."""

    # GIVEN matched normal filters
    normal_filters = BaseSNVFilters.MATCHED_NORMAL_FILTER_NAMES

    # GIVEN all quality SNV filters excluding the matched normal filters
    filter_string_soft_filter_normals = TgaSNVFilters.get_bcftools_filter_string(
        category="quality", analysis_type=AnalysisType.PAIRED, soft_filter_normals=True
    )

    # THEN normal_filters should not exist in the filter string
    for normal_filter in normal_filters:
        assert normal_filter not in filter_string_soft_filter_normals

    # GIVEN all quality SNV filters including the matched normal filters
    filter_string = TgaSNVFilters.get_bcftools_filter_string(
        category="quality", analysis_type=AnalysisType.PAIRED, soft_filter_normals=False
    )

    # THEN normal_filters should exist in the filter string
    for normal_filter in normal_filters:
        assert normal_filter in filter_string


def test_get_bcftools_filter_string():
    """Test bcftools filter string generation."""
    filter_string = TgaSNVFilters.get_bcftools_filter_string(
        category="quality", variant_caller=BioinfoTools.VARDICT
    )
    assert "FILTER~" in filter_string
    assert "Cluster0bp" in filter_string


def test_get_filters_exclude_variantcaller_filters():
    """Test exclusion of variant caller filters."""

    # GIVEN all internal variant caller filters
    internal_variantcaller_filters: List[
        VCFFilter
    ] = BaseSNVFilters.INTERNAL_VARIANT_CALLER_FILTERS
    # GIVEN all quality SNV filters excluding the INTERNAL_VARIANT_CALLER_FILTERS
    retrieved_filters: List[VCFFilter] = TgaSNVFilters.get_filters(
        category="quality", exclude_variantcaller_filters=True
    )

    # WHEN extracting all filter name from both filter lists
    internal_filter_names = {f.filter_name for f in internal_variantcaller_filters}
    retrieved_filter_names = {f.filter_name for f in retrieved_filters}

    # THEN no filter names should overlap in both lists
    assert internal_filter_names.isdisjoint(retrieved_filter_names), (
        "Some filters from internal_variantcaller_filters were found in retrieved_filters: "
        f"{internal_filter_names & retrieved_filter_names}"
    )


# Tests for helper function
def test_get_tag_and_filtername():
    """Test get_tag_and_filtername helper function."""
    # GIVEN 2 dummy VCFFilters
    filters = [
        VCFFilter(tag_value=0.01, filter_name="SWEGENAF", field="INFO"),
        VCFFilter(tag_value=0.1, filter_name="ArtefactFrq", field="INFO"),
    ]
    # WHEN retrieving specific filter tag and filtername
    result = get_tag_and_filtername(filters, "ArtefactFrq")
    # THEN correct tag and filtername should be retrieved
    assert result == [0.1, "ArtefactFrq"]


def test_get_tag_and_filtername_no_match():
    """Test get_tag_and_filtername raises KeyError when no match."""
    # GIVEN dummy VCFFilter list
    filters = [
        VCFFilter(tag_value=0.01, filter_name="SWEGENAF", field="INFO"),
    ]
    # WHEN retrieving tag and filtername from unexisting filtername
    # THEN KeyError exception should be raised
    with pytest.raises(KeyError):
        get_tag_and_filtername(filters, "Nonexistent")


def test_get_tag_and_filtername_multiple_matches():
    """Test get_tag_and_filtername raises ValueError for multiple matches."""
    # GIVEN dummy duplicate VCFFilter names
    filters = [
        VCFFilter(tag_value=0.01, filter_name="Duplicate", field="INFO"),
        VCFFilter(tag_value=0.02, filter_name="Duplicate", field="INFO"),
    ]
    # WHEN retrieving Duplicate filter name
    # THEN ValueError should be raised
    with pytest.raises(ValueError):
        get_tag_and_filtername(filters, "Duplicate")
