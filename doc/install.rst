========
BALSAMIC
========

Bioinformatic Analysis pipeLine for SomAtic MutatIons in Cancer
(**version** = 4.0.0)

.. contents::

Installation
=======

Software Requirments
~~~~~~~~~~~~~~~~~~~~

- Conda: BALSAMIC requires conda 4.3+ for automated installation. For detailed
software and python requirments please see ``requirments.txt`` and
``BALSAMIC/conda/BALSAMIC-base.yaml``

- Singularity: BALSAMIC uses singularity to run vairous parts of the workflow. Either a container has to be built
matching the BALSAMIC version from ``BALSAMIC/containers/BALSAMIC.X.X.X`` or one can pull Singularity container from
Docker using: ``singularity pull docker://hassanf/balsamic path_to_balsamic_container``

- Python 3.6

Manual Installation
~~~~~~~~~~~~~~~~~~~

1. Create a conda environment using ``BALSAMIC/conda/BALSAMIC-base.yaml`` 
2. Install BALSAMIC using ``pip`` within the newly created environment: ``pip install -r requirements.txt -e .``
3. Pull container using Singularity: ``singularity pull docker://hassanf/balsamic path_to_balsamic_container``

Automatic Installation
~~~~~~~~~~~~~~~~~~~~~~

Use ``install.sh`` script, assuming `${CONDA_ENVS_PATH}` is set to the path for conda environment:

::
  
  ./install.sh -s D -v 3.2.2 -p ${CONDA_ENVS_PATH} -c

`-s` set prefix for conda environment name
`-v` change the version tag or branch to append to env name
`-p` path to conda environment location
`-c` use conda 4.6.14 to install packages

In above example, the final conda environment will be named: `D_BALSAMIC-base_3.2.2`

::

  USAGE: ./install.sh [-s _condaprefix -v _balsamic_ver -p _condapath -c]
    1. Conda naming convention: [P,D,S]_[ENVNAME]_[balsamic_version_tag]. P: Production, D: Development, S: Stage
    2. Conda environment prefix: Path to conda env. e.g. /home/user/conda_env/
    
    -s _condaprefix  Conda env name prefix. This will be P or D in the help above
    -v _balsamic_ver Balsamic version tag or branch to install 
    -p _condapath    Conda env path prefix. See point 2 in help above
    -c If set it will use Singularity container instead of local available conda
