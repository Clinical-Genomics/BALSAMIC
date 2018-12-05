========
BALSAMIC
========

Bioinformatic Analysis pipeLine for SomAtic MutatIons in Cancer
(**version** = 2.8.1)

.. contents::

Running BALSAMIC
-----
Prior to run a sample through BALSAMIC, a config file must be generated. In order to create it through an example,
run the following command:

::
    cd test_data
    balsamic config sample \
        --tumor fastq/S1_R_1.fastq.gz \
        --normal fastq/S2_R_1.fastq.gz \
        --panel-bed references/GRCh37/panel/panel.bed \
        --sample-id id1 \
        --analysis-dir ./ \
        --analysis-type paired \
        --reference-config references/reference.json


This will create a directory within test_data directory with the following structure:

