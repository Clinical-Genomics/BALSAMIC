"""Snakemake rules constants."""
from typing import Dict, List

from BALSAMIC.constants.cache import GenomeVersion
from BALSAMIC.constants.analysis import (
    AnalysisType,
    AnalysisWorkflow,
    SequencingType,
    WorkflowSolution,
)

common_cache_rules: List[str] = [
    "snakemake_rules/cache/singularity_containers.rule",
    "snakemake_rules/cache/reference_genome_index.rule",
    "snakemake_rules/cache/reference_download.rule",
]

hg_cache_rules: List[str] = common_cache_rules + [
    "snakemake_rules/cache/cadd.rule",
    "snakemake_rules/cache/delly.rule",
    "snakemake_rules/cache/refseq.rule",
    "snakemake_rules/cache/reference_vcf.rule",
    "snakemake_rules/cache/vep.rule",
]

canfam_cache_rules: List[str] = common_cache_rules + [
    "snakemake_rules/cache/refseq_canfam.rule"
]


SNAKEMAKE_RULES: Dict[str, Dict[str, list]] = {
    "common": {
        "qc": [
            "snakemake_rules/quality_control/fastqc.rule",
            "snakemake_rules/quality_control/multiqc.rule",
            "snakemake_rules/quality_control/qc_metrics.rule",
            "snakemake_rules/quality_control/samtools_qc.rule",
        ],
        "align": [
            "snakemake_rules/align/sentieon_alignment.rule",
            "snakemake_rules/align/bam_compress.rule",
        ],
        "varcall": [
            "snakemake_rules/variant_calling/germline_sv.rule",
            "snakemake_rules/variant_calling/sentieon_quality_filter.rule",
            "snakemake_rules/variant_calling/somatic_sv_quality_filter.rule",
        ],
        "annotate": [
            "snakemake_rules/annotation/somatic_snv_annotation.rule",
            "snakemake_rules/annotation/somatic_sv_annotation.rule",
            "snakemake_rules/annotation/somatic_computations.rule",
            "snakemake_rules/annotation/germline_annotation.rule",
            "snakemake_rules/annotation/varcaller_sv_filter.rule",
            "snakemake_rules/annotation/vcf2cytosure_convert.rule",
            "snakemake_rules/annotation/final_vcf_reheader.rule",
        ],
    },
    "single_targeted": {
        "qc": [
            "snakemake_rules/quality_control/fastp_tga.rule",
            "snakemake_rules/quality_control/picard.rule",
            "snakemake_rules/quality_control/sambamba_depth.rule",
            "snakemake_rules/quality_control/mosdepth.rule",
            "snakemake_rules/umi/concatenation_umi.rule",
            "snakemake_rules/umi/qc_umi.rule",
            "snakemake_rules/umi/mergetype_tumor_umi.rule",
            "snakemake_rules/umi/generate_AF_tables.rule",
        ],
        "align": [
            "snakemake_rules/umi/sentieon_umiextract.rule",
            "snakemake_rules/umi/sentieon_consensuscall.rule",
            "snakemake_rules/align/postprocess_bam.rule",
        ],
        "varcall": [
            "snakemake_rules/variant_calling/germline.rule",
            "snakemake_rules/variant_calling/split_bed.rule",
            "snakemake_rules/variant_calling/somatic_cnv_tumor_only_tga.rule",
            "snakemake_rules/variant_calling/somatic_tumor_only.rule",
            "snakemake_rules/variant_calling/somatic_sv_tumor_only.rule",
            "snakemake_rules/umi/sentieon_varcall_tnscope.rule",
        ],
        "annotate": [
            "snakemake_rules/annotation/rankscore.rule",
            "snakemake_rules/annotation/varcaller_filter_tumor_only.rule",
        ],
    },
    "paired_targeted": {
        "qc": [
            "snakemake_rules/quality_control/fastp_tga.rule",
            "snakemake_rules/quality_control/picard.rule",
            "snakemake_rules/quality_control/sambamba_depth.rule",
            "snakemake_rules/quality_control/mosdepth.rule",
            "snakemake_rules/umi/qc_umi.rule",
            "snakemake_rules/quality_control/somalier.rule",
            "snakemake_rules/umi/concatenation_umi.rule",
            "snakemake_rules/umi/mergetype_tumor_umi.rule",
            "snakemake_rules/umi/mergetype_normal_umi.rule",
            "snakemake_rules/quality_control/contest.rule",
            "snakemake_rules/umi/generate_AF_tables.rule",
        ],
        "align": [
            "snakemake_rules/umi/sentieon_umiextract.rule",
            "snakemake_rules/umi/sentieon_consensuscall.rule",
            "snakemake_rules/align/postprocess_bam.rule",
        ],
        "varcall": [
            "snakemake_rules/variant_calling/germline.rule",
            "snakemake_rules/variant_calling/split_bed.rule",
            "snakemake_rules/variant_calling/somatic_tumor_normal.rule",
            "snakemake_rules/variant_calling/somatic_sv_tumor_normal.rule",
            "snakemake_rules/variant_calling/somatic_cnv_tumor_normal_tga.rule",
            "snakemake_rules/umi/sentieon_varcall_tnscope_tn.rule",
        ],
        "annotate": [
            "snakemake_rules/annotation/rankscore.rule",
            "snakemake_rules/annotation/varcaller_filter_tumor_normal.rule",
            "snakemake_rules/annotation/vcfheader_rename.rule",
        ],
    },
    "single_wgs": {
        "qc": [
            "snakemake_rules/quality_control/fastp_wgs.rule",
            "snakemake_rules/quality_control/sentieon_qc_metrics.rule",
            "snakemake_rules/quality_control/picard_wgs.rule",
            "snakemake_rules/quality_control/report.rule",
        ],
        "varcall": [
            "snakemake_rules/variant_calling/sentieon_germline.rule",
            "snakemake_rules/variant_calling/sentieon_split_snv_sv.rule",
            "snakemake_rules/variant_calling/sentieon_t_varcall.rule",
            "snakemake_rules/variant_calling/somatic_sv_tumor_only.rule",
            "snakemake_rules/dragen_suite/dragen_dna.rule",
        ],
        "annotate": [
            "snakemake_rules/annotation/varcaller_wgs_filter_tumor_only.rule",
        ],
    },
    "paired_wgs": {
        "qc": [
            "snakemake_rules/quality_control/fastp_wgs.rule",
            "snakemake_rules/quality_control/sentieon_qc_metrics.rule",
            "snakemake_rules/quality_control/picard_wgs.rule",
            "snakemake_rules/quality_control/report.rule",
            "snakemake_rules/quality_control/somalier.rule",
        ],
        "varcall": [
            "snakemake_rules/variant_calling/sentieon_germline.rule",
            "snakemake_rules/variant_calling/sentieon_split_snv_sv.rule",
            "snakemake_rules/variant_calling/sentieon_tn_varcall.rule",
            "snakemake_rules/variant_calling/somatic_sv_tumor_normal.rule",
        ],
        "annotate": [
            "snakemake_rules/annotation/varcaller_wgs_filter_tumor_normal.rule",
            "snakemake_rules/annotation/vcfheader_rename.rule",
        ],
    },
    "cache": {
        GenomeVersion.HG19: hg_cache_rules,
        GenomeVersion.HG38: hg_cache_rules,
        GenomeVersion.CanFam3: canfam_cache_rules,
    },
}

DELIVERY_RULES: List[str] = [
    # QC
    "multiqc",
    "collect_custom_qc_metrics",
    "cnv_report",
    # Alignment
    "mergeBam_tumor_umiconsensus",
    "mergeBam_normal_umiconsensus",
    "bam_compress_tumor",
    "bam_compress_normal",
    # Germline
    "vcfheader_rename_germline",
    "vep_annotate_germlineVAR_tumor",
    "vep_annotate_germlineVAR_normal",
    # SNVs
    "bcftools_view_split_variant",
    "bcftools_filter_tnscope_research_tumor_only",
    "bcftools_filter_tnscope_research_tumor_normal",
    "bcftools_filter_tnscope_clinical_tumor_only",
    "bcftools_filter_tnscope_clinical_tumor_normal",
    "vardict_merge",
    "bcftools_filter_vardict_research_tumor_only",
    "bcftools_filter_vardict_research_tumor_normal",
    "bcftools_filter_vardict_clinical_tumor_only",
    "bcftools_filter_vardict_clinical_tumor_normal",
    "sentieon_tnscope_umi",
    "sentieon_tnscope_umi_tn",
    "bcftools_filter_TNscope_umi_research_tumor_only",
    "bcftools_filter_TNscope_umi_research_tumor_normal",
    "bcftools_filter_TNscope_umi_clinical_tumor_only",
    "bcftools_filter_TNscope_umi_clinical_tumor_normal",
    # SVs
    "svdb_merge_tumor_only",
    "svdb_merge_tumor_normal",
    "bcftools_filter_sv_research",
    "bcftools_filter_sv_clinical",
    "tiddit_sv_tumor_only",
    "tiddit_sv_tumor_normal",
    # CNVs
    "delly_cnv_tumor_only",
    "delly_cnv_tumor_normal",
    "ascat_tumor_normal",
    "vcf2cytosure_convert_tumor_only",
    "vcf2cytosure_convert_tumor_normal",
    "cnvkit_segment_CNV_research",
    "cnvkit_call_CNV_research",
    "vcf2cytosure_convert",
    "finalize_gens_outputfiles",
    # TMB
    "tmb_calculation",
]