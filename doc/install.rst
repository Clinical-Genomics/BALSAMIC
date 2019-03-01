========
BALSAMIC
========

Bioinformatic Analysis pipeLine for SomAtic MutatIons in Cancer
(**version** = 2.9.8)

.. contents::

Installation
=======

Software Requirments
~~~~~~~~~~~~~~~~~~~~

BALSAMIC requires conda 4.3+ for automated installation. For detailed
software and python requirments please see ``requirments.txt`` and
``BALSAMIC/install/conda_yaml/*yaml``

BALSAMIC requires multiple conda environments. ``install.sh`` will
automatically install BALSAMIC in the following order:

-  Create a python 3.6 conda environment following the conda environment
   naming convention from
   https://github.com/Clinical-Genomics/development/blob/master/conda/conda_conventions.md
   based on the config file: ``BALSAMIC/conda_yaml/BALSAMIC.yaml``
-  Create conda environments required for BALSAMIC to run properly
-  Install BALSAMIC
-  Install gatk
