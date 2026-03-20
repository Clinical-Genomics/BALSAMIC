
CNV Report Documentation
========================

Overview
--------

The CNV report provides a comprehensive overview of copy number variation (CNV)
signals derived from CNVkit and (optionally) PureCN, together with quality
metrics and visualization.

The report consists of:

- Genome-wide CNV plots (from CNVkit)
- Per-chromosome CNV and VAF plots
- Segment-level table (CNVkit and PureCN)
- Gene-region level table
- Quality control metrics
- Optional PON-based (Panel of Normals) interpretation


Input Data
----------

The report is generated from the following main inputs:

- CNVkit bin-level data (``.cnr``)
- CNVkit segment-level data (``.cns``)
- CNVkit initial segments (``.cns``)
- Germline VCF (for VAF visualization)
- Cancer associated genes list
- Cytoband annotation
- Optional:
  - Panel of Normals (PON) (``.cnn``)
  - PureCN LOH/segment output
  - Purity/ploidy summary



Chromosome Plots
----------------

..  figure:: images/cnvplot_example.png
    :width: 1000px


Each chromosome plot consists of two panels:

Top panel (CNV signal)
^^^^^^^^^^^^^^^^^^^^^^

This panel shows:

- Individual bins (log2 ratios)
- Smoothed log2 signal
- CNV segments (CNVkit and PureCN)
- Optional PON-derived signals
- Highlighted genes

Key features:

**Bin signal**
  Each point represents a genomic bin with a log2 ratio.

**Rolling median**
  A smoothed signal (default: 5-bin window) to reduce noise.

**Segments**
  - CNVkit segments: solid lines
  - PureCN segments: dashed lines

**Color coding**
  - Red: amplification
  - Blue: deletion
  - Black: neutral

**PON spread (if available)**
  Shown as a band representing expected variation across normal samples.

**Highlighted genes**
  Genes of interest are emphasized by:
  - Expanded horizontal space
  - Colored points
  - Gene labels


Bottom panel (VAF / B-allele frequency)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This panel shows:

- Variant allele frequencies (VAF)
- Expected reference lines:
  - 0.5 (heterozygous)
  - 1/3 and 2/3 (imbalance)

If PureCN is available:

- LOH regions are shown as shaded areas

Interpretation:

- Tight clustering around 0.5 → balanced alleles
- Split or shifted clusters → allelic imbalance
- Flattening toward 0 or 1 → possible LOH

Gene Highlighting
-----------------

Genes are highlighted based on:

- Presence of CNV signal (amplification/deletion)
- Presence of LOH (from PureCN)
- Minimum number of targets overlapping gene (corresponding to probes in panel design)

Additionally, see specific logic for cancer associated genes and exome samples below:

**Cancer associated genes:** See :ref:`cancer-gene-lists` for details.

Highlighted if:

- Minimum targets overlapping gene > 2

**Not cancer associated genes:**

Highlighted if:

- CNV or LOH evidence from CNVkit or PureCN
- Minimum targets overlapping gene > 2

Exome-specific behavior
^^^^^^^^^^^^^^^^^^^^^^^

When ``--is-exome`` is enabled:

- Only cancer genes are highlighted
- Non-cancer genes are excluded from highlighting

.. _cancer-gene-lists:

Cancer Gene Lists
-----------------

The report uses a predefined cancer gene list to determine which genes
should be highlighted and prioritized in visualizations and tables.

The gene list is provided as a TSV file with the following columns:

- ``GENE``: gene symbol
- ``SOURCE``: origin of the annotation (e.g. ONCOKB, custom lists)


.. note::

   If there is a gene that you expect to see highlighted in the report but is not included,
   you can request its addition to the cancer gene list.

   Please open an issue in the BALSAMIC repository or contact the maintainers by email.

   Including a short justification (e.g. clinical relevance or panel design) helps ensure
   the request can be evaluated and incorporated quickly.


Selection Rules
^^^^^^^^^^^^^^^

All genes listed in the file are included in the cancer gene set.

- Genes from ONCOKB are included
- Genes from custom or other sources are included
- The list is treated as a curated union of all provided sources

In practice, this means the gene list directly defines which genes are
considered biologically or clinically relevant in the report.

Example
^^^^^^^

Given the following input:

.. code-block:: text

    GENE    OCCURRENCE    SOURCE
    TP53    10            ONCOKB
    EGFR    0             ONCOKB
    MYC     5             ONCOKB
    GENE1   0             cust000

All genes in the file will be included:

- TP53
- EGFR
- MYC
- GENE1

Effect on the Report
^^^^^^^^^^^^^^^^^^^^

The cancer gene set influences several parts of the report:

**Gene highlighting**
  - Cancer genes are prioritized for visualization
  - They may be highlighted even when signal is subtle

**Exome mode**
  - Only cancer genes are considered for highlighting
  - Non-cancer genes are excluded from visualization

**Plot scaling**
  - Highlighted cancer genes are visually expanded in chromosome plots
  - This improves readability of clinically relevant regions

Rationale
^^^^^^^^^

The cancer gene list is intended to represent a curated set of genes
of interest, combining:

- Established cancer genes (e.g. ONCOKB)
- Custom or panel-specific genes

This allows flexibility while ensuring that important regions are not
missed due to strict filtering.


Gene Regions
------------

Gene regions represent aggregated signals across bins belonging to a gene.

Creation steps
^^^^^^^^^^^^^^

1. Bins are grouped per gene
2. Bins are ordered by genomic position
3. Adjacent bins are merged into candidate runs
4. Runs are filtered based on signal strength
5. Small gaps may be bridged if signal is consistent
6. Final regions are collapsed into summary rows

.. note::

   If there is no PON available, the bins of a gene are simply collapsed, and annotated with overlapping CNVkit segment information.

Each region includes:

- Genomic span (start/end)
- Number of targets
- Mean log2 signal

PON-based Interpretation
-----------------------

When a Panel of Normals (PON) is provided, additional metrics are computed.

These include:

- PON mean log2 per region
- PON spread (expected noise)
- Z-score-like signal strength

PON signal classification:

- ``strong``: strong deviation from normal (z-score > 5.0)
- ``borderline``: mild deviation ( 5.0 < z-score > 2.0)
- ``noise``: no significant deviation (z-score =< 2.0)

PON indication:

- ``GAIN``: likely amplification ((mean_log2 - mean_log2_pon) > 0.07)
- ``LOSS``: likely deletion ((mean_log2 - mean_log2_pon) < -0.07)


.. note::

   Only "strong" signals are can be classified as GAIN / LOSS, otherwise, gene-region is set to NEUTRAL

   Only genes with a minimum of 8 targets will be considered for this PON based GAIN / LOSS indication


Important:

- PON interpretation is only available if a PON file is provided


Segment Table
-------------

The segment table combines CNVkit and PureCN results.

Columns include:

- Chromosome and genomic coordinates
- Segment size
- CNVkit log2 and copy number
- PureCN copy number and LOH annotations
- Cytoband
- Genes in segment (see note below!)

CNV calls are standardized:

- Amplification
- Deletion
- Neutral


.. note::

   Only genes with a minimum of 2 targets are shown in the segment table

   For **exome** this is further limited to only show cancer-associated genes


Sex-aware interpretation
^^^^^^^^^^^^^^^^^^^^^^^^

Copy number interpretation accounts for sex chromosomes:

- X and Y are interpreted differently depending on sample sex
- Prevents misclassification of normal sex chromosome states


Gene Region Table
-----------------

The gene-region table summarizes CNV signals at gene level.

Includes:

- Gene name
- Number of targets
- Mean log2
- CNVkit and PureCN calls
- PON-based metrics (if available)

This table is useful for:

- Seeing indications for focal amplifications/deletions that may be missed by official callers

.. warning::

   The gene-region analysis is based on a custom aggregation method and has not
   been formally validated for clinical use.

   The signals presented here should be interpreted as **supportive evidence only**
   and must not be used as standalone or definitive CNV calls.

   Any findings from this table should be confirmed using established CNV callers
   (e.g. CNVkit, PureCN) and/or orthogonal methods before being considered for
   clinical interpretation.


Quality Metrics
---------------

The report includes summary QC metrics:

Log2 noise
^^^^^^^^^^

- Derived from adjacent bin differences
- Similar to CNVkit DLRSpread
- Lower values indicate cleaner signal

Filtered targets
^^^^^^^^^^^^^^^^

- Number of bins removed relative to PON
- Indicates how much of the original panel design could be used for CNV analysis


Limitations
-----------

- CNV detection depends on coverage quality and tumor purity
- Small focal events may be missed in low-coverage regions
- PON interpretation depends on quality and composition of normal samples
- LOH detection requires reliable PureCN results

