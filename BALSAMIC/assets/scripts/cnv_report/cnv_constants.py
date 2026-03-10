from dataclasses import dataclass, field
from typing import Mapping, Sequence

CHR = "chr"


@dataclass(frozen=True)
class TableSpec:
    column_order: Sequence[str]
    float_columns: Sequence[str]
    decimals: int = 3
    renames: Mapping[str, str] | None = None
    sort_keys: tuple[str, str, str] = ("chr", "region_start", "region_end")
    descriptions: Mapping[str, str] = field(default_factory=dict)


RENAME_COLS = {
    "chr": "Chr",
    "start": "Start",
    "end": "End",
    "region_start": "Start",
    "region_end": "End",
    "pon_region_z": "PON deviation score",
    "pon_region_indication": "PON-based indication",
    "pon_region_signal": "PON deviation signal",
    "cnvkit_cnv_call": "CNVkit segment call",
    "purecn_cnv_call": "PureCN segment call",
    "purecn_type": "PureCN segment LOH-type",
    "gene.symbol": "Gene",
    "baf_maf": "BAF | MAF",
}


SEGMENT_COLUMN_DESCRIPTIONS = {
    # --- Genomic coordinates ---
    "chr": "Chromosome identifier.",
    "start": "Start coordinate (bp) of the CN segment.",
    "end": "End coordinate (bp) of the CN segment.",
    "cytoband": "Cytogenetic band annotation for the gene region.",
    "segment_size": "Size of the segment: [end - start]",
    "caller": "The caller which determined the segment, CNVkit / PureCN",
    # --- Gene identity ---
    "gene.symbol": "HGNC gene symbols overlapping segment (minimum 2 probe targets required for a gene to be listed)",
    # --- Bin aggregation metrics ---
    "n.targets": "Number of target bins overlapping with CN segment.",
    "log2": "log2 copy-ratio from CNVkit or PureCN (not adjusted for purity or ploidy).",
    "baf_maf": "Mean B-allele frequency (BAF) for the CNVkit segment / Observed minor allele frequency (MAF) in the PureCN segment.",
    "cnv_call": "CNV call from PureCN or CNVkit (AMPLIFICATION / DELETION / NEUTRAL) determined based on total copy-numbers predicted from each tool.",
    "purecn_type": "PureCN call annotation string for the segment.",
    # --- CNVkit segment details ---
    "cnvkit_seg_depth": "Mean sequencing read depth across all targets belonging to the CNVkit segment.",
    "cnvkit_adjusted_log2": "Segment-level log2 copy-ratio from CNVkit, adjusted for purity and ploidy.",
    "cnvkit_seg_cn": "Total copy number estimate from CNVkit.",
    "cnvkit_seg_cn1": "Estimated minor allele copy number from CNVkit.",
    "cnvkit_seg_cn2": "Estimated major allele copy number from CNVkit.",
    # --- PureCN segment details --
    "purecn_num_snps": "Number of SNPs used by PureCN in the segment.",
    "purecn_C": "Total copy number estimate from PureCN.",
    "purecn_M": "Minor allele copy number estimate from PureCN.",
    "purecn_M_flagged": "True if PureCN flagged the minor allele estimate as unreliable or special case.",
}

GENE_REGION_COLUMN_DESCRIPTIONS = {
    # --- Genomic coordinates ---
    "chr": "Chromosome identifier (without 'chr' prefix).",
    "region_start": "Start coordinate (bp) of the aggregated gene region.",
    "region_end": "End coordinate (bp) of the aggregated gene region.",
    "cytoband": "Cytogenetic band annotation for the gene region.",
    # --- Gene identity ---
    "gene.symbol": "HGNC gene symbol.",
    "is_cancer_gene": "True if the gene is included in the configured cancer gene set.",
    # --- Bin aggregation metrics ---
    "n.targets": "Number of target bins contributing to this gene-level region.",
    "mean_log2": "Mean log2 copy-ratio across bins overlapping the gene region.",
    "min_log2": "Minimum log2 copy-ratio observed among bins in the gene region.",
    "max_log2": "Maximum log2 copy-ratio observed among bins in the gene region.",
    "depth_mean": "Mean sequencing depth across bins in the gene region.",
    # --- CNVkit gene-level call ---
    "cnvkit_cnv_call": "Gene-level CNV classification derived from CNVkit segmentation (e.g., DELETION, AMPLIFICATION, NEUTRAL).",
    # --- PureCN gene-level call ---
    "purecn_cnv_call": "Gene-level CNV classification derived from PureCN.",
    "purecn_type": "PureCN call annotation string for the overlapping segment; non-empty indicates LOH-type classification.",
    # --- CNVkit segment details ---
    "cnvkit_seg_start": "Start coordinate (bp) of the overlapping CNVkit segment.",
    "cnvkit_seg_end": "End coordinate (bp) of the overlapping CNVkit segment.",
    "cnvkit_seg_log2": "Segment-level log2 copy-ratio from CNVkit.",
    "cnvkit_seg_raw_log2": "Raw (pre-adjustment) segment-level log2 copy-ratio from CNVkit.",
    "cnvkit_seg_baf": "Mean B-allele frequency (BAF) for the overlapping CNVkit segment.",
    "cnvkit_seg_cn": "Total copy number estimate from CNVkit.",
    "cnvkit_seg_cn1": "Estimated minor allele copy number from CNVkit.",
    "cnvkit_seg_cn2": "Estimated major allele copy number from CNVkit.",
    # --- PureCN segment details ---
    "purecn_seg_start": "Start coordinate (bp) of the overlapping PureCN segment.",
    "purecn_seg_end": "End coordinate (bp) of the overlapping PureCN segment.",
    "purecn_seg_mean_log2": "Segment-level log2 copy-ratio from PureCN.",
    "purecn_num_snps": "Number of SNPs used by PureCN in the segment.",
    "purecn_maf_observed": "Observed minor allele frequency (MAF) in the PureCN segment.",
    "purecn_C": "Total copy number estimate from PureCN.",
    "purecn_M": "Minor allele copy number estimate from PureCN.",
    "purecn_M_flagged": "True if PureCN flagged the minor allele estimate as unreliable or special case.",
    # --- Exon overlap ---
    "exons_overlapping_cnvkit_segment": "Number of exons overlapping the CNVkit segment for this gene.",
    "exons_overlapping_gene_region": "Number of exons overlapping the aggregated gene region.",
    # --- PON bin/region-level statistics ---
    "pon_mean_log2": "Mean PON log2 value across bins in the gene-region.",
    "pon_mean_spread": "Mean PON spread across bins in the gene-region.",
    "pon_region_effect": "Difference between sample log2 and PON mean at gene region level.",
    "pon_region_z": "Z-score of gene region-level deviation relative to PON variability.",
    "pon_region_signal": "Signal classification of region-level PON deviation (noisy / borderline / strong).",
    "pon_region_indication": "Interpretation of gene regionp-level deviation based on PON (e.g., GAIN/LOSS/NEUTRAL).",
}


SEGMENT_TABLE_SPEC = TableSpec(
    column_order=[
        "chr",
        "start",
        "end",
        "cytoband",
        "cnv_call",
        "purecn_type",
        "log2",
        "baf_maf",
        "n.targets",
        "segment_size",
        "cnvkit_adjusted_log2",
        "cnvkit_seg_cn",
        "cnvkit_seg_cn1",
        "cnvkit_seg_cn2",
        "cnvkit_seg_depth",
        "purecn_C",
        "purecn_M",
        "purecn_M_flagged",
        "purecn_num_snps",
        "caller",
        "gene.symbol",
    ],
    float_columns=[
        "log2",
        "baf_maf",
        "cnvkit_adjusted_log2",
        "cnvkit_seg_cn",
        "cnvkit_seg_cn1",
        "cnvkit_seg_cn2",
        "cnvkit_seg_depth",
        "purecn_C",
        "purecn_M",
    ],
    decimals=3,
    sort_keys=("chr", "start", "end"),
    renames=RENAME_COLS,
    descriptions=SEGMENT_COLUMN_DESCRIPTIONS,
)

GENE_TABLE_SPEC = TableSpec(
    column_order=[
        "chr",
        "region_start",
        "region_end",
        "cytoband",
        "gene.symbol",
        "n.targets",
        "mean_log2",
        "min_log2",
        "max_log2",
        "cnvkit_cnv_call",
        "purecn_cnv_call",
        "purecn_type",
        "cnvkit_seg_start",
        "cnvkit_seg_end",
        "cnvkit_seg_log2",
        "cnvkit_seg_raw_log2",
        "cnvkit_seg_baf",
        "purecn_seg_start",
        "purecn_seg_end",
        "purecn_seg_mean_log2",
        "purecn_num_snps",
        "purecn_maf_observed",
        "cnvkit_seg_cn",
        "cnvkit_seg_cn1",
        "cnvkit_seg_cn2",
        "purecn_C",
        "purecn_M",
        "purecn_M_flagged",
        "is_cancer_gene",
        "depth_mean",
        "pon_mean_log2",
        "pon_mean_spread",
        "pon_region_effect",
        "pon_region_z",
        "pon_region_signal",
        "pon_region_indication",
    ],
    float_columns=[
        "mean_log2",
        "min_log2",
        "max_log2",
        "depth_mean",
        "pon_mean_log2",
        "pon_mean_spread",
        "pon_region_effect",
        "pon_region_z",
        "cnvkit_seg_log2",
        "cnvkit_seg_raw_log2",
        "cnvkit_seg_baf",
        "purecn_seg_mean_log2",
        "purecn_maf_observed",
    ],
    decimals=3,
    renames=RENAME_COLS,
    descriptions=GENE_REGION_COLUMN_DESCRIPTIONS,
)


PURECN_WARNING_TEXT = (
    "PureCN failed to identify a reliable final purity/ploidy model for this sample. "
    "The displayed purity and ploidy values are default fallback values and should not "
    "be interpreted as inferred sample-specific estimates. "
    "Furthermore, CNVkit Copy Number determination which is based on PureCN purity is also unreliable. "
    "Finally; PureCN CNV segments will not be produced, and LOH will not be indicated. "
)
