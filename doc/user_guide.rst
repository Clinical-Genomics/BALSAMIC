========
BALSAMIC
========

Bioinformatic Analysis pipeLine for SomAtic MutatIons in Cancer
(**version** = 3.1.1)

.. contents::

Running BALSAMIC
-----
Prior to run a sample through BALSAMIC, a config file must be generated. In order to create it through an example,
run the following command:

::

    cd test_data
    balsamic config case \
        --tumor fastq/S1_R_1.fastq.gz \
        --normal fastq/S2_R_1.fastq.gz \
        --panel-bed references/GRCh37/panel/panel.bed \
        --case-id id1 \
        --analysis-dir ./ \
        --analysis-type paired \
        --output-config id1_analysis.json \
        --reference-config references/reference.json


This will create a directory within test_data directory with the following structure:

::

    id1
      ├── BALSAMIC_run
      ├── analysis
      ├── fastq
      │   ├── S1_R_1.fastq.gz -> ../fastq/S1_R_1.fastq.gz
      │   ├── S1_R_2.fastq.gz -> ../fastq/S1_R_2.fastq.gz
      │   ├── S2_R_1.fastq.gz -> ../fastq/S2_R_1.fastq.gz
      │   └── S2_R_2.fastq.gz -> ../fastq/S2_R_2.fastq.gz
      ├── id1_analysis.json
      ├── id1_analysis.json_BALSAMIC_3.1.1_graph.pdf
      ├── logs
      └── scripts


And to test run balsamic:

::

  balsamic run -s id1/id1_analysis.json

BALSAMIC is always in dry-run mode. If you want to run an actual config file, make sure you add `-f` flag.
