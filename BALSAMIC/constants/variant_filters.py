from BALSAMIC.models.params import VCFFilter
from BALSAMIC.constants.analysis import (
    AnalysisType,
    BioinfoTools,
)
from typing import List, Optional, Literal, Set
from enum import Enum


class BaseSNVFilters:
    INTERNAL_VARIANT_CALLER_FILTERS = [
        VCFFilter(
            filter_name="Cluster0bp",
            Description="Two variants are within 0 bp",
            variant_caller=BioinfoTools.VARDICT,
        ),
        VCFFilter(
            filter_name="InGap",
            Description="The variant is in the deletion gap, thus likely false positive",
            variant_caller=BioinfoTools.VARDICT,
        ),
        VCFFilter(
            filter_name="InIns",
            Description="The variant is adjacent to an insertion variant",
            variant_caller=BioinfoTools.VARDICT,
        ),
        VCFFilter(
            filter_name="MSI12",
            Description="Variant in MSI region with 12 non-monomer MSI or 13 monomer MSI",
            variant_caller=BioinfoTools.VARDICT,
        ),
        VCFFilter(
            filter_name="NM4.5",
            Description="Mean mismatches in reads >= 4.5, thus likely false positive",
            variant_caller=BioinfoTools.VARDICT,
        ),
        VCFFilter(
            filter_name="SN1.5",
            Description="Signal to Noise Less than 1.5",
            variant_caller=BioinfoTools.VARDICT,
        ),
        VCFFilter(
            filter_name="f0.001",
            Description="Allele frequency < 0.001",
            variant_caller=BioinfoTools.VARDICT,
        ),
        VCFFilter(
            filter_name="p8",
            Description="Mean Position in Reads Less than 8",
            variant_caller=BioinfoTools.VARDICT,
        ),
        VCFFilter(
            filter_name="pSTD",
            Description="Position in Reads has STD of 0",
            variant_caller=BioinfoTools.VARDICT,
        ),
        VCFFilter(
            filter_name="q22.5",
            Description="Mean Base Quality Below 22.5",
            variant_caller=BioinfoTools.VARDICT,
        ),
        VCFFilter(
            filter_name="Bias",
            Description="Strand Bias",
            variant_caller=BioinfoTools.VARDICT,
        ),
        VCFFilter(
            filter_name="AMPBIAS",
            Description="Indicate the variant has amplicon bias",
            variant_caller=BioinfoTools.VARDICT,
            analysis_type=AnalysisType.SINGLE,
        ),
        VCFFilter(
            filter_name="LongMSI",
            Description="The somatic variant is flanked by long A/T (>=14)",
            variant_caller=BioinfoTools.VARDICT,
            analysis_type=AnalysisType.SINGLE,
        ),
        VCFFilter(
            filter_name="Q10",
            Description="Mean Mapping Quality Below 10",
            variant_caller=BioinfoTools.VARDICT,
            analysis_type=AnalysisType.SINGLE,
        ),
        VCFFilter(
            filter_name="d3",
            Description="Allele frequency < 0.001",
            variant_caller=BioinfoTools.VARDICT,
            analysis_type=AnalysisType.SINGLE,
        ),
        VCFFilter(
            filter_name="v2",
            Description="Var Depth < 2",
            variant_caller=BioinfoTools.VARDICT,
            analysis_type=AnalysisType.SINGLE,
        ),
        VCFFilter(
            filter_name="DIFF0.2",
            Description="Non-somatic or LOH and allele frequency difference < 0.2",
            variant_caller=BioinfoTools.VARDICT,
            analysis_type=AnalysisType.PAIRED,
        ),
        VCFFilter(
            filter_name="LongAT",
            Description="The somatic variant is flanked by long A/T (>=14)",
            variant_caller=BioinfoTools.VARDICT,
            analysis_type=AnalysisType.PAIRED,
        ),
        VCFFilter(
            filter_name="P0.9",
            Description="Not significant with p-value > 0.9",
            variant_caller=BioinfoTools.VARDICT,
            analysis_type=AnalysisType.PAIRED,
        ),
        VCFFilter(
            filter_name="d5",
            Description="Total Depth < 5",
            variant_caller=BioinfoTools.VARDICT,
            analysis_type=AnalysisType.PAIRED,
        ),
        VCFFilter(
            filter_name="v3",
            Description="Var Depth < 3",
            variant_caller=BioinfoTools.VARDICT,
            analysis_type=AnalysisType.PAIRED,
        ),
        VCFFilter(
            filter_name="Q0",
            Description="Mean Mapping Quality Below 0",
            variant_caller=BioinfoTools.VARDICT,
            analysis_type=AnalysisType.PAIRED,
        ),
        VCFFilter(
            filter_name="t_lod_fstar",
            Description="Tumor does not meet likelihood threshold",
            variant_caller=BioinfoTools.TNSCOPE,
        ),
        VCFFilter(
            filter_name="low_t_alt_frac",
            Description="Site filtered due to low alt allele fraction",
            variant_caller=BioinfoTools.TNSCOPE,
        ),
        VCFFilter(
            filter_name="germline_risk",
            Description="Evidence indicates this site is germline, not somatic",
            variant_caller=BioinfoTools.TNSCOPE,
            analysis_type=AnalysisType.PAIRED,
        ),
    ]

    MATCHED_NORMAL_FILTER_NAMES: List[str] = [
        "germline_risk",
        "in_normal",
    ]
    DEFAULT_SOFT_FILTERS: List[str] = [
        "HighOccurrenceFrq",
    ]

    @classmethod
    def filter_criteria(
        cls,
        category: Literal["clinical", "research", "quality"],
        analysis_type: Optional[Enum] = None,
        variant_caller: Optional[Enum] = None,
        exclude_variantcaller_filters: Optional[bool] = False,
        exome: Optional[bool] = False,
    ) -> List[VCFFilter]:
        """
        Shared filtering logic to get filters based on criteria.

        Args:
            category (Literal["clinical", "research", "quality"]): The filter category to use.
            analysis_type (Optional[Enum]): Filter based on analysis type (default: None).
            variant_caller (Optional[Enum]): Filter based on variant caller (default: None).
            exclude_variantcaller_filters (Optional[bool]): If True, excludes the variantcaller filters.
            exome (Optional[bool]): Filter based on exome sequencing (default: False).
        Returns:
            List[VCFFilter]: A list of matching filter objects.
        """
        # Get the filters for the category
        filters = getattr(cls, category)

        # Helper function to check if a filter matches the criteria
        def filter_matches(f: VCFFilter) -> bool:
            return (
                (
                    analysis_type is None
                    or getattr(f, "analysis_type", None) in {None, analysis_type}
                )
                and (
                    variant_caller is None
                    or getattr(f, "variant_caller", None) in {None, variant_caller}
                )
                and (exome is None or getattr(f, "exome", None) in {None, exome})
            )

        # Filter the filters based on the matching function
        matching_filters = [f for f in filters if filter_matches(f)]

        # If category is "quality", include VARIANT_CALLER_HARDFILTERS
        if category == "quality" and not exclude_variantcaller_filters:
            matching_filters += [
                f for f in cls.INTERNAL_VARIANT_CALLER_FILTERS if filter_matches(f)
            ]

        return matching_filters

    @classmethod
    def get_bcftools_filter_string(
        cls,
        category: Literal["clinical", "research", "quality"],
        analysis_type: Optional[Enum] = None,
        variant_caller: Optional[Enum] = None,
        soft_filter_normals: Optional[bool] = None,
        exome: Optional[bool] = False,
    ) -> str:
        """
        Get a set of filter names based on various attributes.

        Args:
            category (Literal["clinical", "research", "quality"]): The filter category to use.
            analysis_type (Optional[Enum]): Filter based on analysis type (default: None).
            variant_caller (Optional[Enum]): Filter based on variant caller (default: None).
            soft_filter_normals (Optional[bool]): If True, removes filters in MATCHED_NORMAL_FILTER_NAMES.
            exome (Optional[bool]): Filter based on exome sequencing (default: False).
        Returns:
            str: bcftools filter string
        """
        # Use the shared filtering logic and extract filter names
        filters = cls.filter_criteria(
            category=category,
            analysis_type=analysis_type,
            variant_caller=variant_caller,
            exome=exome,
        )

        # Exclude default soft filters
        filters = [f for f in filters if f.filter_name not in cls.DEFAULT_SOFT_FILTERS]

        # Exclude filters if soft_filter_normals is set
        if soft_filter_normals:
            filters = [
                f
                for f in filters
                if f.filter_name not in cls.MATCHED_NORMAL_FILTER_NAMES
            ]

        # Extract filter_names
        filter_names = [f.filter_name for f in filters]

        # Format as BCFTools-compatible filter string
        return " || ".join(
            [f'FILTER~"{filter_name}"' for filter_name in sorted(filter_names)]
        )

    @classmethod
    def get_filters(
        cls,
        category: Literal["clinical", "research", "quality"],
        analysis_type: Optional[Enum] = None,
        variant_caller: Optional[Enum] = None,
        exclude_variantcaller_filters: Optional[bool] = True,
        exome: Optional[bool] = False,
    ) -> List[VCFFilter]:
        """
        Get a list of filters matching the specified attributes.

        Args:
            category (Literal["clinical", "research", "quality"]): The filter category to use.
            analysis_type (Optional[Enum]): Filter based on analysis type (default: None).
            variant_caller (Optional[Enum]): Filter based on variant caller (default: None).
            exclude_variantcaller_filters (Optional[bool]): If True, excludes the variantcaller filters.
            exome (Optional[bool]): Filter based on exome sequencing (default: False).
        Returns:
            List[VCFFilter]: A list of matching filter objects.
        """
        # Use the shared filtering logic and return filter objects
        return cls.filter_criteria(
            category,
            analysis_type,
            variant_caller,
            exclude_variantcaller_filters=exclude_variantcaller_filters,
            exome=exome,
        )


class WgsSNVFilters(BaseSNVFilters):
    research = [
        VCFFilter(tag_value=0.01, filter_name="SWEGENAF", field="INFO"),
        VCFFilter(tag_value=0.001, filter_name="balsamic_high_pop_freq", field="INFO"),
    ]
    clinical = research + [
        VCFFilter(tag_value=0.007, filter_name="Frq", field="INFO"),
        VCFFilter(tag_value=0.1, filter_name="ArtefactFrq", field="INFO"),
    ]
    quality = [
        VCFFilter(
            tag_value=0.3,
            filter_name="in_normal",
            field="FORMAT",
            analysis_type=AnalysisType.PAIRED,
        ),
        VCFFilter(tag_value=3, filter_name="balsamic_low_tumor_ad", field="FORMAT"),
        VCFFilter(tag_value=10, filter_name="balsamic_low_tumor_dp", field="FORMAT"),
        VCFFilter(tag_value=0.05, filter_name="balsamic_low_af", field="FORMAT"),
        VCFFilter(
            tag_value=20,
            filter_name="balsamic_low_quality_scores",
            field="FORMAT",
            analysis_type=AnalysisType.SINGLE,
        ),
        VCFFilter(
            tag_value=0,
            filter_name="balsamic_low_strand_read_counts",
            field="FORMAT",
            analysis_type=AnalysisType.SINGLE,
        ),
        VCFFilter(
            tag_value=3,
            filter_name="balsamic_high_strand_oddsratio",
            field="INFO",
            analysis_type=AnalysisType.SINGLE,
        ),
        VCFFilter(
            tag_value=4,
            filter_name="balsamic_high_strand_oddsratio",
            field="INFO",
            analysis_type=AnalysisType.PAIRED,
        ),
    ]


class TgaSNVFilters(BaseSNVFilters):
    research = [
        VCFFilter(
            filter_name="MERGED", Description="SNV Merged with neighboring variants"
        ),
        VCFFilter(tag_value=0.01, filter_name="SWEGENAF", field="INFO"),
        VCFFilter(tag_value=0.005, filter_name="balsamic_high_pop_freq", field="INFO"),
    ]
    clinical = research + [
        VCFFilter(tag_value=0.01, filter_name="Frq", field="INFO"),
        VCFFilter(tag_value=0.1, filter_name="ArtefactFrq", field="INFO"),
        VCFFilter(tag_value=0.3, filter_name="HighOccurrenceFrq", field="INFO"),
    ]
    quality = [
        VCFFilter(
            tag_value=0.3,
            filter_name="in_normal",
            field="FORMAT",
            analysis_type=AnalysisType.PAIRED,
        ),
        VCFFilter(tag_value=0.005, filter_name="balsamic_low_af", field="INFO"),
        VCFFilter(tag_value=5, filter_name="balsamic_low_tumor_ad", field="INFO"),
        VCFFilter(
            tag_value=50,
            filter_name="balsamic_low_tumor_dp",
            field="INFO",
            exome=False,
        ),
        VCFFilter(
            tag_value=20,
            filter_name="balsamic_low_tumor_dp",
            field="INFO",
            exome=True,
        ),
        VCFFilter(
            tag_value=30,
            filter_name="balsamic_low_mq",
            field="INFO",
            variant_caller=BioinfoTools.VARDICT,
        ),
        VCFFilter(
            tag_value=20,
            filter_name="balsamic_low_quality_scores",
            field="INFO",
            variant_caller=BioinfoTools.TNSCOPE,
            analysis_type=AnalysisType.SINGLE,
        ),
        VCFFilter(
            tag_value=2.7,
            filter_name="balsamic_high_strand_oddsratio",
            field="INFO",
            variant_caller=BioinfoTools.TNSCOPE,
            analysis_type=AnalysisType.SINGLE,
        ),
        VCFFilter(
            tag_value=3,
            filter_name="balsamic_high_strand_oddsratio",
            field="INFO",
            variant_caller=BioinfoTools.TNSCOPE,
            analysis_type=AnalysisType.PAIRED,
        ),
        VCFFilter(
            tag_value=12,
            filter_name="balsamic_high_tnscope_rpa",
            field="INFO",
            variant_caller=BioinfoTools.TNSCOPE,
        ),
    ]


class TgaUmiSNVFilters(BaseSNVFilters):
    research = [
        VCFFilter(tag_value=0.01, filter_name="SWEGENAF", field="INFO"),
        VCFFilter(tag_value=0.02, filter_name="balsamic_high_pop_freq", field="INFO"),
    ]
    clinical = research + [
        VCFFilter(tag_value=0.01, filter_name="Frq", field="INFO"),
        VCFFilter(tag_value=0.1, filter_name="ArtefactFrq", field="INFO"),
    ]
    quality = [
        VCFFilter(
            tag_value=0.3,
            filter_name="in_normal",
            field="FORMAT",
            analysis_type=AnalysisType.PAIRED,
        ),
    ]


def get_tag_and_filtername(filters: List[VCFFilter], filter_name: str) -> List[str]:
    """
    Retrieve the filter name and tag value for a specific filter.

    Args:
        filters (List[VCFFilter]): A list of filter objects to search within.
        filter_name (str): The name of the filter to match.

    Returns:
        List[str]: A list containing the tag value and the filter name.

    Raises:
        ValueError: If more than one filter matches the filter name.
        KeyError: If no filter matches the filter name.
    """
    matched_filters = [f for f in filters if f.filter_name == filter_name]

    if len(matched_filters) > 1:
        raise ValueError(
            f"Multiple filters found with the filter name '{filter_name}'. "
        )
    elif not matched_filters:
        raise KeyError(f"No filter found with the filter name '{filter_name}'.")

    # Return the filter name and tag value for the matched filter
    return [matched_filters[0].tag_value, matched_filters[0].filter_name]


# Manta bcftools filters
MANTA_FILTER_SETTINGS = {
    "low_pr_sr_count": {
        "tag_value": 4,
        "filter_name": "low_pr_sr_count",
        "field": "FORMAT",
    },
    "varcaller_name": "Manta",
    "filter_type": "general",
    "analysis_type": "tumor_only,tumor_normal",
    "description": "Bcftools filters to set frequency and minimum read support for SV calls",
}

# Configuration for SVDB settings:
SVDB_FILTER_SETTINGS = {
    "swegen_sv_freq": {
        "tag_value": 0.02,
        "filter_name": "SWEGENAF",
        "field": "INFO",
    },
    "loqusdb_clinical_sv_freq": {
        "tag_value": 0.02,
        "filter_name": "Frq",
        "field": "INFO",
    },
    "varcaller_name": "svdb",
    "filter_type": "general",
    "analysis_type": "tumor_only,tumor_normal",
    "description": "General purpose filters used for filtering svdb merged sv",
}
