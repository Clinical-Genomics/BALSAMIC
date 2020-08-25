============
Installation
============

This section describes steps to install BALSAMIC (**version** = 4.1.0)

.. contents::

Software Requirments
~~~~~~~~~~~~~~~~~~~~

- Conda >=version 4.5.0: For detailed software and python requirments please see ``requirments.txt`` and ``BALSAMIC/conda/balsamic.yaml``

- Singularity >=version 3.0.0: BALSAMIC uses singularity to run vairous parts of the workflow. Either a container has to
  be built matching the BALSAMIC version from ``BALSAMIC/containers/Dockerfile.latest`` or one can pull Singularity
container from Docker using: ``singularity pull path_container_file docker://hassanf/balsamic:tag`` 

- Python 3.6

Manual Installation
~~~~~~~~~~~~~~~~~~~

1. Create a conda environment using ``BALSAMIC/conda/balsamic.yaml``:

::

conda env create --file BALSAMIC/conda/balsamic.yaml 


You can also set conda environment prefix via: `--prefix`
 
2. Install BALSAMIC using ``pip`` within the newly created environment: ``pip install -r requirements.txt -e .``

3. Pull container using Singularity: ``singularity pull path_container_file docker://hassanf/balsamic``


Example:
If you'd like to install release 5.0.0 the instruction will look like below:

::

# Create a conda env: balsamic_base
conda env create --file BALSAMIC/conda/balsamic.yaml --quiet --force --prefx balsamic_base

# Activate conda environment
source activate balsamic_base

# Pull container for release_v5.0.0
singularity pull balsamic_release_v5.0.0 docker://hassanf/balsamic:release_v5.0.0


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
