from dataclasses import dataclass, field
from typing import Mapping, Sequence

CHR = "chr"
GENE = "gene.symbol"
TARGETS = "n.targets"


@dataclass(frozen=True)
class GeneRegionConfig:
    # ------------------------------------------------------------------
    # Initial gene-region detection
    # ------------------------------------------------------------------

    # Minimum number of bins required in a gene before attempting to
    # detect internal gene regions. Genes with few targets are skipped because
    # runs cannot be detected robustly.
    min_gene_targets: int = 8

    # Minimum number of consecutive bins required to form a candidate run.
    # Prevents short noisy stretches from becoming gene regions.
    min_run_bins: int = 4

    # Minimum absolute per-bin z-score required for a bin to participate
    # in a candidate run:
    #
    #     z = (log2 - pon_log2) / pon_spread
    #
    # Bins with |z| below this threshold break a run.
    z_bin_thresh: float = 1.5

    # Minimum aggregated run-level z-score required to accept a candidate
    # run as a provisional gene region. This summarizes the entire run
    # using its mean deviation and expected noise.
    z_run_thresh: float = 3.0

    # Median smoothing window applied to bin-level z-scores before run
    # detection. Helps suppress isolated noisy bins that would otherwise
    # split runs.
    smooth_window: int = 3

    # ------------------------------------------------------------------
    # Region cleanup / merging
    # ------------------------------------------------------------------

    # Maximum number of bins allowed in a "bridge" segment between two
    # larger regions with similar effects. Used to merge patterns like:
    #
    #     strong run → small neutral gap → strong run
    #
    # when the gap is short enough.
    max_bridge_bins: int = 4

    # Maximum difference in mean effect between two runs for them to be
    # considered similar enough to merge across a bridge segment.
    #
    # Effect is defined as:
    #
    #     effect = log2 - pon_log2
    bridge_delta: float = 0.12

    # ------------------------------------------------------------------
    # Region scoring and CNV indication
    # ------------------------------------------------------------------

    # Minimum number of bins required for a region to be eligible for
    # PON-based statistical scoring and CNV indication.
    #
    # Regions smaller than this are treated as unreliable and will not
    # produce gain/loss calls.
    min_region_targets_for_call: int = 5

    # Lower threshold for classifying the PON z-score signal strength.
    #
    # z < noise_lt       → "noise"
    # noise_lt ≤ z < borderline_lt → "borderline"
    # z ≥ borderline_lt  → "strong"
    pon_signal_noise_lt: float = 2.0

    # Upper threshold separating "borderline" from "strong" signals.
    pon_signal_borderline_lt: float = 5.0

    # Log2 deviation threshold above which a region is interpreted as
    # a gain when the signal is considered strong.
    pon_gain_gt: float = 0.07

    # Log2 deviation threshold below which a region is interpreted as
    # a loss when the signal is considered strong.
    pon_loss_lt: float = -0.07

    # Neutral call string
    pon_neutral_call = "NEUTRAL"

    # Gain call string
    pon_gain_call = "GAIN"

    # Loss call string
    pon_loss_call = "LOSS"


@dataclass(frozen=True)
class ChromosomePlotConfig:
    backbone_factor: float
    neutral_target_factor: float
    min_gene_targets: int
    min_gene_targets_cancer: int
    highlight_only_cancer: bool
    log2_rolling_window: int = 5
    base_label_offset: float = 1.5
    y_abs_max: float = 3.0
    key: tuple[str, str, str] = ("chr", "start", "end")


PANEL_PLOT_CONFIG = ChromosomePlotConfig(
    backbone_factor=0.4,
    neutral_target_factor=0.4,
    min_gene_targets=2,
    min_gene_targets_cancer=2,
    highlight_only_cancer=False,
)

EXOME_PLOT_CONFIG = ChromosomePlotConfig(
    backbone_factor=0.07,
    neutral_target_factor=0.07,
    min_gene_targets=4,  # this will not be used as highlight cancer is True
    min_gene_targets_cancer=4,
    highlight_only_cancer=True,
)


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
    "purecn_type": "PureCN LOH-type",
    GENE: "Gene",
    "baf_maf": "BAF | MAF",
    "segment_size": "Segment size (KB)",
}

pnp = "(adjusted for purity and ploidy)"
not_pnp = "(NOT adjusted for purity or ploidy)"

SEGMENT_COLUMN_DESCRIPTIONS = {
    # --- Genomic coordinates ---
    "chr": "Chromosome identifier.",
    "start": "Start coordinate (bp) of the CN segment.",
    "end": "End coordinate (bp) of the CN segment.",
    "cytoband": "Cytogenetic band annotation for the gene region.",
    "segment_size": "Size of the segment in Kbp.",
    "caller": "The caller which determined the segment, CNVkit / PureCN",
    # --- Gene identity ---
    GENE: "HGNC gene symbols overlapping segment (minimum 2 probe targets required for a gene to be listed)",
    # --- Bin aggregation metrics ---
    TARGETS: "Number of target bins overlapping with CN segment.",
    "log2": f"log2 copy-ratio from CNVkit or PureCN {not_pnp}.",
    "baf_maf": "Mean B-allele frequency (BAF) for the CNVkit segment | Median observed minor allele frequency (MAF) of heterozygous SNPs in the PureCN segment.",
    "cnv_call": f"CNV call from PureCN or CNVkit (AMPLIFICATION / DELETION / NEUTRAL) determined based on total copy-numbers predicted from each tool.",
    "purecn_type": "PureCN LOH type if M = 0 (such as WHOLE ARM COPY-NEUTRAL LOH, COPY-NEUTRAL LOH)",
    # --- CNVkit segment details ---
    "cnvkit_seg_depth": "Mean sequencing read depth across all targets belonging to the CNVkit segment.",
    "cnvkit_adjusted_log2": f"Segment-level log2 copy-ratio from CNVkit {pnp}.",
    "cnvkit_seg_cn": f"Total copy number estimate from CNVkit {pnp}.",
    "cnvkit_seg_cn1": f"Estimated minor allele copy number from CNVkit {pnp}.",
    "cnvkit_seg_cn2": f"Estimated major allele copy number from CNVkit {pnp}.",
    # --- PureCN segment details --
    "purecn_num_snps": "Number of SNPs in this segment informative for LOH detection from PureCN.",
    "purecn_C": f"Total copy number for segment from PureCN {pnp}.",
    "purecn_M": f"Minor allele copy number from PureCN {pnp}.",
    "purecn_M_flagged": "True if PureCN flagged the minor allele estimate as unreliable (due to small number of variants or ambiguous call)",
}

GENE_REGION_COLUMN_DESCRIPTIONS = {
    # --- Genomic coordinates ---
    "chr": "Chromosome identifier (without 'chr' prefix).",
    "region_start": "Start coordinate (bp) of the aggregated gene region.",
    "region_end": "End coordinate (bp) of the aggregated gene region.",
    "cytoband": "Cytogenetic band annotation for the gene region.",
    # --- Gene identity ---
    GENE: "HGNC gene symbol.",
    "is_cancer_gene": "True if the gene is included in the configured cancer gene set.",
    # --- Bin aggregation metrics ---
    TARGETS: "Number of target bins in this gene region. Genes are split into multiple rows when log2 ratios relative to the PON suggest distinct within-gene copy number segments (GAIN or LOSS).",
    "mean_log2": f"Mean log2 copy number ratio across bins overlapping the gene region {not_pnp}.",
    "min_log2": f"Minimum log2 copy number ratio observed among bins in the gene region {not_pnp}.",
    "max_log2": f"Maximum log2 copy number ratio observed among bins in the gene region {not_pnp}.",
    "depth_mean": "Mean sequencing depth across bins in the gene region.",
    # --- CNVkit gene-level call ---
    "cnvkit_cnv_call": "CNV call from overlapping CNVkit segment (e.g., DELETION, AMPLIFICATION, NEUTRAL).",
    # --- PureCN gene-level call ---
    "purecn_cnv_call": "CNV call from overlapping PureCN segment (e.g., DELETION, AMPLIFICATION, NEUTRAL).",
    "purecn_type": "PureCN LOH type if M = 0 (such as WHOLE ARM COPY-NEUTRAL LOH, COPY-NEUTRAL LOH).",
    # --- CNVkit segment details ---
    "cnvkit_seg_start": "Start coordinate (bp) of the overlapping CNVkit segment.",
    "cnvkit_seg_end": "End coordinate (bp) of the overlapping CNVkit segment.",
    "cnvkit_adjusted_log2": f"Segment-level log2 copy-ratio from CNVkit {pnp}.",
    "cnvkit_seg_raw_log2": f"Raw segment-level log2 copy-ratio from CNVkit {not_pnp}.",
    "cnvkit_seg_baf": "Mean B-allele frequency (BAF) for the overlapping CNVkit segment.",
    "cnvkit_seg_cn": f"Total copy number estimate for the overlapping CNVkit segment {pnp}.",
    "cnvkit_seg_cn1": f"Estimated minor allele copy number for the overlapping CNVkit segment {pnp}.",
    "cnvkit_seg_cn2": f"Estimated major allele copy number for the overlapping CNVkit segment {pnp}.",
    # --- PureCN segment details ---
    "purecn_seg_start": "Start coordinate (bp) of the overlapping PureCN segment.",
    "purecn_seg_end": "End coordinate (bp) of the overlapping PureCN segment.",
    "purecn_seg_mean_log2": f"Segment mean of copy number log2-ratios from PureCN {not_pnp}.",
    "purecn_num_snps": "Number of SNPs in this segment informative for LOH detection for the overlapping PureCN segment",
    "purecn_maf_observed": "Median observed minor allele frequency (MAF) of heterozygous SNPs for the overlapping PureCN segment.",
    "purecn_C": f"Total copy number for segment for the overlapping PureCN segment {pnp}.",
    "purecn_M": f"Minor allele copy number for the overlapping PureCN segment {pnp}.",
    "purecn_M_flagged": "For the overlapping PureCN segment: True if PureCN flagged the minor allele estimate as unreliable (due to small number of variants or ambiguous call).",
    # --- PON bin/region-level statistics ---
    "pon_mean_log2": "Mean PON log2 value across bins in the gene-region.",
    "pon_mean_spread": "Mean PON spread across bins in the gene-region.",
    "pon_region_effect": "Difference between sample log2 and PON mean at sub-gene region level.",
    "pon_region_z": "Deviation score relative to PON for sub-gene region-level (z-score like).",
    "pon_region_signal": "Signal classification for sub-gene region-level PON deviation (noisy / borderline / strong).",
    "pon_region_indication": "PON-based GAIN/LOSS/NEUTRAL indication for sub-gene region.",
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
        TARGETS,
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
        GENE,
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
        GENE,
        TARGETS,
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
