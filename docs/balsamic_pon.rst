Panel of Normals (PON)
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

    balsamic config pon --case-id <CASE_ID> --balsamic-cache </path/reference_cache/> --analysis-dir </path/analysis/> --fastq-path </path/fastq/> --panel-bed </path/panel.bed>

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

