========
BALSAMIC
========

Bioinformatic Analysis pipeLine for SomAtic MutatIons in Cancer
(**version** = 3.1.1)

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
-  ``install.sh`` has following options:

::

  USAGE: ./install.sh [-s _condaprefix -d _condadate -p _condapath -c]
    1. Conda naming convention: [P,D,S]_[ENVNAME]_%DATE. P: Production, D: Development
    2. Conda environment prefix: Path to conda env. e.g. /home/user/conda_env/
    
    -s _condaprefix  Conda env name prefix. This will be P, D, or S in the help above. 
    -d _condadate    Conda env name suffix. This will be a suffix, by default it will be current date: yymmdd 
    -p _condapath    Conda env path prefix. See point 2 in help above.
    -c If set it will use Singularity container for conda instead

-  Example command to install BALSAMIC and its environments:

::

  ./install.sh -s D -d 190905 -p path_to_conda_env -c

