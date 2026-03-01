CHR = "chr"


from dataclasses import dataclass, field
from typing import Mapping, Sequence


@dataclass(frozen=True)
class TableSpec:
    column_order: Sequence[str]
    float_columns: Sequence[str]
    decimals: int = 3
    renames: Mapping[str, str] | None = None
    sort_keys: tuple[str, str, str] = ("chr", "region_start", "region_end")
    descriptions: Mapping[str, str] = field(default_factory=dict)


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
        "exons_overlapping_cnvkit_segment",
        "is_cancer_gene",
        "depth_mean",
        "pon_gene_mean_log2",
        "pon_gene_mean_spread",
        "pon_gene_effect",
        "pon_gene_z",
        "pon_gene_direction",
        "pon_gene_significance",
        "pon_gene_indication",
        "pon_mean_log2",
        "pon_mean_spread",
        "pon_chunk_effect",
        "pon_chunk_z",
        "pon_chunk_significance",
        "pon_chunk_indication",
        "exons_overlapping_gene_region",
    ],
    float_columns=[
        "seg_baf",
        "mean_log2",
        "min_log2",
        "max_log2",
        "pon_gene_mean_log2",
        "pon_gene_mean_spread",
        "pon_gene_effect",
        "pon_gene_z",
        "depth_mean",
        "pon_mean_log2",
        "pon_mean_spread",
        "pon_chunk_effect",
        "pon_chunk_z",
        "cnvkit_seg_log2",
        "cnvkit_seg_raw_log2",
        "cnvkit_seg_baf",
        "purecn_seg_mean_log2",
        "purecn_maf_observed",
    ],
    decimals=3,
    renames={"n_targets": "n.targets"},
    descriptions={
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
        "purecn_seg_mean": "Segment-level log2 copy-ratio from PureCN.",
        "purecn_num_snps": "Number of SNPs used by PureCN in the segment.",
        "purecn_maf_observed": "Observed minor allele frequency (MAF) in the PureCN segment.",
        "purecn_C": "Total copy number estimate from PureCN.",
        "purecn_M": "Minor allele copy number estimate from PureCN.",
        "purecn_M_flagged": "True if PureCN flagged the minor allele estimate as unreliable or special case.",
        # --- Exon overlap ---
        "exons_overlapping_cnvkit_segment": "Number of exons overlapping the CNVkit segment for this gene.",
        "exons_overlapping_gene_region": "Number of exons overlapping the aggregated gene region.",
        # --- PON gene-level statistics ---
        "pon_gene_mean_log2": "Mean PON log2 value across bins overlapping the gene region.",
        "pon_gene_mean_spread": "Mean PON spread (variability) across bins overlapping the gene region.",
        "pon_gene_effect": "Difference between sample log2 and PON mean at gene level.",
        "pon_gene_z": "Z-score of the gene-level effect relative to PON variability.",
        "pon_gene_direction": "Direction of deviation relative to PON (e.g., GAIN or LOSS).",
        "pon_gene_significance": "Significance classification of the gene-level PON deviation.",
        "pon_gene_indication": "Interpretation of gene-level deviation based on PON (e.g., GAIN/LOSS/NEUTRAL).",
        # --- PON bin/chunk-level statistics ---
        "pon_mean_log2": "Mean PON log2 value across bins in the chunk or region.",
        "pon_mean_spread": "Mean PON spread across bins in the chunk or region.",
        "pon_chunk_effect": "Difference between sample log2 and PON mean at chunk level.",
        "pon_chunk_z": "Z-score of chunk-level deviation relative to PON variability.",
        "pon_chunk_significance": "Significance classification of chunk-level PON deviation.",
        "pon_chunk_indication": "Interpretation of chunk-level deviation based on PON (e.g., GAIN/LOSS/NEUTRAL).",
    },
)
