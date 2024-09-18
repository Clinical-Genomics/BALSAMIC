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
        "misc": ["snakemake_rules/misc/sleep.rule"],
        "qc": [
            "snakemake_rules/quality_control/fastqc.rule",
            "snakemake_rules/quality_control/multiqc.rule",
            "snakemake_rules/quality_control/qc_metrics.rule",
            "snakemake_rules/quality_control/picard_common.rule",
            "snakemake_rules/quality_control/sentieon_qc_metrics.rule",
        ],
        "report": [
            "snakemake_rules/report/generate_pdf.rule",
            "snakemake_rules/report/merge_pdfs.rule",
        ],
        "align": [
            "snakemake_rules/align/bam_compress.rule",
        ],
        "varcall": [
            "snakemake_rules/variant_calling/snv_quality_filter.rule",
            "snakemake_rules/variant_calling/tnscope_post_process.rule",
        ],
        "annotate": [
            "snakemake_rules/annotation/somatic_snv_annotation.rule",
            "snakemake_rules/annotation/somatic_sv_annotation.rule",
            "snakemake_rules/annotation/somatic_computations.rule",
            "snakemake_rules/annotation/germline_annotation.rule",
            "snakemake_rules/annotation/varcaller_sv_filter.rule",
            "snakemake_rules/annotation/vcf2cytosure_convert.rule",
            "snakemake_rules/annotation/final_vcf_reheader.rule",
            "snakemake_rules/annotation/rankscore.rule",
        ],
    },
    "single_targeted": {
        "qc": [
            "snakemake_rules/quality_control/fastp_tga.rule",
            "snakemake_rules/quality_control/picard.rule",
            "snakemake_rules/quality_control/sambamba_depth.rule",
            "snakemake_rules/quality_control/mosdepth.rule",
            "snakemake_rules/concatenation/concatenation.rule",
            "snakemake_rules/umi/qc_umi.rule",
            "snakemake_rules/umi/generate_AF_tables.rule",
            "snakemake_rules/quality_control/samtools_qc_tga.rule",
        ],
        "align": [
            "snakemake_rules/align/tga_sentieon_alignment.rule",
            "snakemake_rules/align/tga_bam_postprocess.rule",
            "snakemake_rules/umi/umi_sentieon_alignment.rule",
            "snakemake_rules/umi/sentieon_consensuscall.rule",
        ],
        "varcall": [
            "snakemake_rules/variant_calling/extend_bed.rule",
            "snakemake_rules/variant_calling/cnvkit_preprocess.rule",
            "snakemake_rules/variant_calling/germline_tga.rule",
            "snakemake_rules/variant_calling/somatic_cnv_tumor_only_tga.rule",
            "snakemake_rules/variant_calling/somatic_sv_tumor_only_tga.rule",
            "snakemake_rules/umi/sentieon_varcall_tnscope.rule",
            "snakemake_rules/umi/modify_tnscope_infofield_umi",
            "snakemake_rules/variant_calling/snv_t_varcall_tga.rule",
            "snakemake_rules/variant_calling/somatic_sv_postprocess_and_filter_tumor_only.rule",
            "snakemake_rules/variant_calling/merge_snv_vcfs.rule",
            "snakemake_rules/variant_calling/vardict_pre_and_postprocessing.rule",
        ],
        "annotate": [
            "snakemake_rules/annotation/varcaller_filter_tumor_only.rule",
            "snakemake_rules/annotation/varcaller_filter_tumor_only_umi.rule",
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
            "snakemake_rules/concatenation/concatenation.rule",
            "snakemake_rules/umi/generate_AF_tables.rule",
            "snakemake_rules/quality_control/samtools_qc_tga.rule",
        ],
        "align": [
            "snakemake_rules/align/tga_sentieon_alignment.rule",
            "snakemake_rules/align/tga_bam_postprocess.rule",
            "snakemake_rules/umi/sentieon_consensuscall.rule",
            "snakemake_rules/umi/umi_sentieon_alignment.rule",
        ],
        "varcall": [
            "snakemake_rules/variant_calling/extend_bed.rule",
            "snakemake_rules/variant_calling/cnvkit_preprocess.rule",
            "snakemake_rules/variant_calling/germline_tga.rule",
            "snakemake_rules/variant_calling/somatic_sv_tumor_normal_tga.rule",
            "snakemake_rules/variant_calling/somatic_cnv_tumor_normal_tga.rule",
            "snakemake_rules/umi/sentieon_varcall_tnscope_tn.rule",
            "snakemake_rules/umi/modify_tnscope_infofield_umi",
            "snakemake_rules/variant_calling/snv_tn_varcall_tga.rule",
            "snakemake_rules/variant_calling/somatic_sv_postprocess_and_filter_tumor_normal.rule",
            "snakemake_rules/variant_calling/merge_snv_vcfs.rule",
            "snakemake_rules/variant_calling/vardict_pre_and_postprocessing.rule",
        ],
        "annotate": [
            "snakemake_rules/annotation/varcaller_filter_tumor_normal.rule",
            "snakemake_rules/annotation/varcaller_filter_tumor_normal_umi.rule",
            "snakemake_rules/annotation/vcfheader_rename.rule",
            "snakemake_rules/annotation/msi_tumor_normal.rule",
        ],
    },
    "single_wgs": {
        "qc": [
            "snakemake_rules/quality_control/fastp_wgs.rule",
            "snakemake_rules/quality_control/picard_wgs.rule",
            "snakemake_rules/quality_control/samtools_qc_wgs.rule",
        ],
        "align": [
            "snakemake_rules/align/wgs_bam_postprocess.rule",
            "snakemake_rules/align/wgs_sentieon_alignment.rule",
        ],
        "varcall": [
            "snakemake_rules/variant_calling/germline_wgs.rule",
            "snakemake_rules/variant_calling/sentieon_t_varcall_wgs.rule",
            "snakemake_rules/variant_calling/somatic_sv_tumor_only_wgs.rule",
            "snakemake_rules/dragen_suite/dragen_dna.rule",
            "snakemake_rules/variant_calling/somatic_sv_postprocess_and_filter_tumor_only.rule",
        ],
        "annotate": [
            "snakemake_rules/annotation/varcaller_filter_tumor_only_wgs.rule",
        ],
    },
    "paired_wgs": {
        "qc": [
            "snakemake_rules/quality_control/fastp_wgs.rule",
            "snakemake_rules/quality_control/picard_wgs.rule",
            "snakemake_rules/quality_control/somalier.rule",
            "snakemake_rules/quality_control/samtools_qc_wgs.rule",
        ],
        "align": [
            "snakemake_rules/align/wgs_bam_postprocess.rule",
            "snakemake_rules/align/wgs_sentieon_alignment.rule",
        ],
        "varcall": [
            "snakemake_rules/variant_calling/germline_wgs.rule",
            "snakemake_rules/variant_calling/sentieon_tn_varcall_wgs.rule",
            "snakemake_rules/variant_calling/somatic_sv_tumor_normal_wgs.rule",
            "snakemake_rules/variant_calling/somatic_sv_postprocess_and_filter_tumor_normal.rule",
        ],
        "annotate": [
            "snakemake_rules/annotation/varcaller_filter_tumor_normal_wgs.rule",
            "snakemake_rules/annotation/vcfheader_rename.rule",
            "snakemake_rules/annotation/msi_tumor_normal.rule",
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
    "modify_tnscope_infofield",
    "modify_tnscope_infofield_umi",
    "gatk_update_vcf_sequence_dictionary",
    "bcftools_filter_tnscope_research_tumor_only",
    "bcftools_filter_tnscope_research_tumor_normal",
    "bcftools_filter_tnscope_clinical_tumor_only",
    "bcftools_filter_tnscope_clinical_tumor_normal",
    "sentieon_tnscope_umi",
    "sentieon_tnscope_umi_tn",
    "bcftools_filter_TNscope_umi_research_tumor_only",
    "bcftools_filter_TNscope_umi_research_tumor_normal",
    "bcftools_filter_TNscope_umi_clinical_tumor_only",
    "bcftools_filter_TNscope_umi_clinical_tumor_normal",
    "genmod_score_snvs",
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
    "cnvpytor_tumor_only",
    "vcf2cytosure_convert_tumor_only",
    "vcf2cytosure_convert_tumor_normal",
    "cnvkit_segment_CNV_research",
    "cnvkit_call_CNV_research",
    "vcf2cytosure_convert",
    "finalize_gens_outputfiles",
    # TMB
    "tmb_calculation",
    # CNV report
    "merge_cnv_pdf_reports",
]
