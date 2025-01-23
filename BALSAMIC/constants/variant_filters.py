
from BALSAMIC.models.params import VCFFilter
from BALSAMIC.constants.analysis import SequencingType, AnalysisWorkflow, AnalysisType, BioinfoTools
from typing import List, Optional, Literal
from enum import Enum

class BaseSNVFilters:
    VARIANT_CALLER_HARDFILTERS = [VCFFilter(filter_name="Cluster0bp", Description="Two variants are within 0 bp", variant_caller=BioinfoTools.VARDICT),
        VCFFilter(filter_name="InGap", Description="The variant is in the deletion gap, thus likely false positive", variant_caller=BioinfoTools.VARDICT),
        VCFFilter(filter_name="InIns", Description="The variant is adjacent to an insertion variant", variant_caller=BioinfoTools.VARDICT),
        VCFFilter(filter_name="MSI12", Description="Variant in MSI region with 12 non-monomer MSI or 13 monomer MSI", variant_caller=BioinfoTools.VARDICT),
        VCFFilter(filter_name="NM4.5", Description="Mean mismatches in reads >= 4.5, thus likely false positive", variant_caller=BioinfoTools.VARDICT),
        VCFFilter(filter_name="SN1.5", Description="Signal to Noise Less than 1.5", variant_caller=BioinfoTools.VARDICT),
        VCFFilter(filter_name="f0.001", Description="Allele frequency < 0.001", variant_caller=BioinfoTools.VARDICT),
        VCFFilter(filter_name="p8", Description="Mean Position in Reads Less than 8", variant_caller=BioinfoTools.VARDICT),
        VCFFilter(filter_name="pSTD", Description="Position in Reads has STD of 0", variant_caller=BioinfoTools.VARDICT),
        VCFFilter(filter_name="AMPBIAS", Description="Indicate the variant has amplicon bias", variant_caller=BioinfoTools.VARDICT,analysis_type=AnalysisType.SINGLE),
        VCFFilter(filter_name="LongMSI", Description="The somatic variant is flanked by long A/T (>=14)", variant_caller=BioinfoTools.VARDICT, analysis_type=AnalysisType.SINGLE),
        VCFFilter(filter_name="Q10", Description="Mean Mapping Quality Below 10", variant_caller=BioinfoTools.VARDICT, analysis_type=AnalysisType.SINGLE),
        VCFFilter(filter_name="d3", Description="Allele frequency < 0.001", variant_caller=BioinfoTools.VARDICT,analysis_type=AnalysisType.SINGLE),
        VCFFilter(filter_name="v2", Description="Var Depth < 2", variant_caller=BioinfoTools.VARDICT, analysis_type=AnalysisType.SINGLE),
        VCFFilter(filter_name="Bias", Description="Strand Bias", variant_caller=BioinfoTools.VARDICT, analysis_type=AnalysisType.PAIRED),
        VCFFilter(filter_name="DIFF0.2", Description="Non-somatic or LOH and allele frequency difference < 0.2", variant_caller=BioinfoTools.VARDICT, analysis_type=AnalysisType.PAIRED),
        VCFFilter(ilter_name="LongAT", Description="The somatic variant is flanked by long A/T (>=14)", variant_caller=BioinfoTools.VARDICT, analysis_type=AnalysisType.PAIRED),
        VCFFilter(filter_name="MAF0.05", Description="Matched sample has AF > 0.05, thus not somatic", variant_caller=BioinfoTools.VARDICT, analysis_type=AnalysisType.PAIRED),
        VCFFilter(filter_name="P0.9",Description="Not significant with p-value > 0.9", variant_caller=BioinfoTools.VARDICT, analysis_type=AnalysisType.PAIRED),
        VCFFilter(filter_name="d5", Description="Total Depth < 5", variant_caller=BioinfoTools.VARDICT, analysis_type=AnalysisType.PAIRED),
        VCFFilter(filter_name="v3", Description="Var Depth < 3", variant_caller=BioinfoTools.VARDICT, analysis_type=AnalysisType.PAIRED),
        VCFFilter(filter_name="t_lod_fstar", Description="Tumor does not meet likelihood threshold", variant_caller=BioinfoTools.TNSCOPE),
        VCFFilter(filter_name="low_t_alt_frac", Description="Site filtered due to low alt allele fraction", variant_caller=BioinfoTools.TNSCOPE),
        VCFFilter(filter_name="germline_risk", Description="Evidence indicates this site is germline, not somatic", variant_caller=BioinfoTools.TNSCOPE, analysis_type=AnalysisType.PAIRED),
        ]

    MATCHED_NORMAL_FILTER_NAMES: List[str] = [
        "germline_risk",
        "high_normal_tumor_af_frac",
        "MAF0.05",
    ]

    @classmethod
    def get_filter_names(
        cls,
        category: Literal["clinical", "research", "quality"],
        analysis_workflow: Optional[Enum] = None,
        analysis_type: Optional[Enum] = None,
        sequencing_type: Optional[Enum] = None,
        variant_caller: Optional[Enum] = None,
        soft_filter_normals: Optional[bool] = None,
    ) -> List[str]:
        """
        Get a list of filter names based on various attributes.

        Args:
            category (Literal["clinical", "research", "quality"]): The filter category to use.
            analysis_workflow (Optional[Enum]): Filter based on analysis workflow (default: None).
            analysis_type (Optional[Enum]): Filter based on analysis type (default: None).
            sequencing_type (Optional[Enum]): Filter based on sequencing type (default: None).
            variant_caller (Optional[Enum]): Filter based on variant caller (default: None).
            soft_filter_normals (Optional[bool]): If True, removes filters in MATCHED_NORMAL_FILTER_NAMES.

        Returns:
            List[str]: A list of filter names.
        """
        filters = getattr(cls, category)
        filtered_names = [
            f.filter_name
            for f in filters
            if (analysis_workflow is None or getattr(f, "analysis_workflow", None) == analysis_workflow)
               and (analysis_type is None or getattr(f, "analysis_type", None) == analysis_type)
               and (sequencing_type is None or getattr(f, "sequencing_type", None) == sequencing_type)
               and (variant_caller is None or getattr(f, "variant_caller", None) == variant_caller)
        ]

        # If category is "quality", also include VARIANT_CALLER_HARDFILTERS
        if category == "quality":
            filtered_names += [
                f.filter_name
                for f in cls.VARIANT_CALLER_HARDFILTERS
                if (analysis_workflow is None or getattr(f, "analysis_workflow", None) == analysis_workflow)
                   and (analysis_type is None or getattr(f, "analysis_type", None) == analysis_type)
                   and (sequencing_type is None or getattr(f, "sequencing_type", None) == sequencing_type)
                   and (variant_caller is None or getattr(f, "variant_caller", None) == variant_caller)
            ]

        if soft_filter_normals:
            filtered_names = [filter_name for filter_name in filtered_names if filter_name not in cls.MATCHED_NORMAL_FILTER_NAMES]

        return filtered_names
class WGS_SNV_Filters(BaseSNVFilters):
    clinical = [
        VCFFilter(tag_value=0.01, filter_name="Frq", field="INFO"),
        VCFFilter(tag_value=0.1, filter_name="ArtefactFrq", field="INFO"),
    ]
    research = [
        VCFFilter(tag_value=0.01, filter_name="SWEGENAF", field="INFO"),
        VCFFilter(tag_value=0.001, filter_name="balsamic_high_pop_freq", field="INFO"),
    ]
    quality = [
        VCFFilter(tag_value=0.3, filter_name="high_normal_tumor_af_frac", field="FORMAT", analysis_type=AnalysisType.PAIRED),
        VCFFilter(tag_value=3, filter_name="balsamic_low_tumor_ad", field="FORMAT"),
        VCFFilter(tag_value=10, filter_name="balsamic_low_tumor_dp", field="FORMAT"),
        VCFFilter(tag_value=0.05, filter_name="balsamic_low_af", field="FORMAT"),
        VCFFilter(tag_value=20, filter_name="balsamic_low_quality_scores", field="FORMAT", analysis_type=AnalysisType.SINGLE),
        VCFFilter(tag_value=0, filter_name="balsamic_low_strand_read_counts", field="FORMAT", analysis_type=AnalysisType.SINGLE),
        VCFFilter(tag_value=3, filter_name="balsamic_high_strand_oddsratio", field="INFO", analysis_type=AnalysisType.SINGLE),
    ]


class TGA_SNV_Filters(BaseSNVFilters):
    clinical = [
        VCFFilter(tag_value=0.01, filter_name="Frq", field="INFO"),
        VCFFilter(tag_value=0.1, filter_name="ArtefactFrq", field="INFO"),
    ]
    research = [
        VCFFilter(tag_value=0.01, filter_name="SWEGENAF", field="INFO"),
        VCFFilter(tag_value=0.005, filter_name="balsamic_high_pop_freq", field="INFO"),
        VCFFilter(tag_value=0.02, filter_name="balsamic_high_pop_freq", field="INFO", analysis_workflow=AnalysisWorkflow.BALSAMIC_UMI),
    ]
    quality = [
        VCFFilter(tag_value=0.3, filter_name="high_normal_tumor_af_frac", field="FORMAT", analysis_type=AnalysisType.PAIRED),
        VCFFilter(tag_value=0.005, filter_name="balsamic_low_af", field="INFO"),
        VCFFilter(tag_value=5, filter_name="balsamic_low_tumor_ad", field="INFO"),
        VCFFilter(tag_value=50, filter_name="balsamic_low_tumor_dp", field="INFO"),
        VCFFilter(tag_value=20, filter_name="balsamic_low_tumor_dp", field="INFO", sequencing_type=SequencingType.WES),
        VCFFilter(tag_value=30, filter_name="balsamic_low_mq", field="INFO", variant_caller=BioinfoTools.VARDICT),
        VCFFilter(tag_value=20, filter_name="balsamic_low_quality_scores", field="INFO", variant_caller=BioinfoTools.TNSCOPE, analysis_type=AnalysisType.SINGLE),
        VCFFilter(tag_value=2.7, filter_name="balsamic_high_strand_oddsratio", field="INFO", variant_caller=BioinfoTools.TNSCOPE,  analysis_type=AnalysisType.SINGLE),
    ]


class TGA_UMI_SNV_Filters(BaseSNVFilters):
    clinical = [
        VCFFilter(tag_value=0.01, filter_name="Frq", field="INFO"),
        VCFFilter(tag_value=0.1, filter_name="ArtefactFrq", field="INFO"),
    ]
    research = [
        VCFFilter(tag_value=0.01, filter_name="SWEGENAF", field="INFO"),
        VCFFilter(tag_value=0.02, filter_name="balsamic_high_pop_freq", field="INFO"),
    ]
    quality = [
        VCFFilter(tag_value=0.3, filter_name="high_normal_tumor_af_frac", field="FORMAT", analysis_type=AnalysisType.PAIRED),
    ]



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
