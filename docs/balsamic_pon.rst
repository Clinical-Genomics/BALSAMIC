Panel of Normals (PON)
======================

Currently two PON-methods are implemented in BALSAMIC to correct for biases and normalise coverage values:

- For producing more accurate CNV variant-calls using ``CNVkit`` for TGA cases.

- To produce normalised CN-profiles for WGS cases visualised in ``GENS``.


CNVkit PON
======================

BALSAMIC provides a functionality to generate a Panel of Normals (PON) for more accurate copy-number filtering of false positives and that can be used as an input for the ``CNVkit`` variant caller.

For a more detailed PON use case, please refer to the following documentation:

- `CNVkit`_
- `Illumina DRAGEN`_
- `GATK`_

.. _CNVkit: https://cnvkit.readthedocs.io/en/stable/pipeline.html#paired-or-pooled-normals
.. _Illumina DRAGEN: https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/GPipelineVarCalNorm_fDG.htm
.. _GATK: https://gatk.broadinstitute.org/hc/en-us/articles/360035890631-Panel-of-Normals-PON-

PON Generation
--------------

When creating a new PON reference file, the next steps have to be followed:

1. Identify the samples to be included in the PON and add their ``fastq`` files to the ``fastq`` directory

.. note::

    One needs to fetch normal samples coming from the same origin, tissue or blood

2. Generate the ``<CASE_ID>_PON.json`` configuration file:

.. code-block::

    balsamic config pon --pon-workflow CNVkit --case-id <CASE_ID> --balsamic-cache </path/reference_cache/> --analysis-dir </path/analysis/> --fastq-path </path/fastq/> --panel-bed </path/panel.bed>

3. Run the BALSAMIC PON workflow:

.. code-block::

    balsamic run analysis -s </path/analysis/<CASE_ID>/<CASE_ID>_PON.json -r


4. Check for the PON reference finish and output files:

.. code-block::

    /path/analysis/analysis_PON_finish
    /path/analysis/cnv/<panel_name>_CNVkit_PON_reference_<version>.cnn

Using the PON during analysis
-----------------------------

BALSAMIC can use a PON reference file if its provided while running CNVkit analysis:

.. code-block::

    balsamic config case --case-id <CASE_ID> --pon-cnn /path/analysis/cnv/<panel_name>_CNVkit_PON_reference_<version>.cnn --balsamic-cache </path/reference_cache/> --analysis-dir </path/analysis/> --panel-bed </path/panel.bed> --tumor-path </path/tumor.fastq>


.. note::

    In the absence of a PON reference file, CNVkit is capable of generating a flat reference (tumor-only) or normal reference (tumor-normal) file on its own to correct for GC content and regional coverage

GENS PON
======================

In order to produce an accurate CN-profile to visualise in GENS you need to create 2 PONs one for each gender (see instructions below).

The original instructions for how to create this PON, and which has been implemented in this BALSAMIC workflow can be found on the Clinical-Genomics-Lund GENS-repository:

- `Clinical-Genomics-Lund-GENS`_

.. _Clinical-Genomics-Lund-GENS: https://github.com/Clinical-Genomics-Lund/gens

To create the PON using the GENS PON creation workflow you can follow the guide below.

PON Generation
--------------

To create a GENS PON using the BALSAMIC workflow you need to follow these steps:

1. Create a genome-interval file.

**Note:**

    These are the genome bins within which the coverage will be calculated, and consequently is the lowest resolution of viewing the CN-profile.

This is the setting we used:

.. code-block::

    gatk PreprocessIntervals --reference [ref] --bin-length 100 --interval-merging-rule OVERLAPPING_ONLY -O human_g1k_v37_gens_targets_preprocessed_100bp.interval_list


2. Identify the samples to be included in the PON and add or link their ``fastq`` files to the ``fastq`` directory

**Note:**

    It is recommended to include approximately 100 samples of the same gender, using the same library preparation and sequencing method as your intended analysis-samples.

2. Generate the ``<CASE_ID>_PON.json`` configuration file:

.. code-block::

    balsamic config pon --pon-creation-type <[GENS_female,GENS_male]> --genome-interval <[path-to-file-from-step1]> --case-id <CASE_ID> --balsamic-cache </path/reference_cache/> --analysis-dir </path/analysis/> --fastq-path </path/fastq/> --panel-bed </path/panel.bed>

3. Run the BALSAMIC PON workflow:

**Note:**
    If you are following these instructions using 100 WGS samples, you require access to compute-nodes with a lot of memory (one of our jobs crashed at 117GB).

.. code-block::

    balsamic run analysis -s </path/analysis/<CASE_ID>/<CASE_ID>_PON.json -r

This workflow runs trimming and alignment for all samples to be included in the PON. Calculates coverages in bins using ``GATK CollectReadCounts`` then creates the PON using all read-counts with the tool ``GATK CreateReadCountPanelOfNormals``.

4. Check for the PON output files:

.. code-block::

    /path/analysis/analysis_PON_finish
    /path/analysis/cnv/gens_pon_100bp.<GENDER>.<VERSION>.hdf5

Using the PON during analysis
-----------------------------

This PON is a required input in order to produce the final output-files to be loaded into the GENS platform.

How to run a case using this PON and to activate GENS for your WGS analysis you are referred to this page:

`Using GENS for WGS <https://balsamic.readthedocs.io/en/latest/balsamic_sv_cnv.html>`_.