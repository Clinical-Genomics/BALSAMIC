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

    conda env create --file BALSAMIC/conda/balsamic.yaml -n name_of_environment 


You can also set conda environment prefix via ``--prefix`` option. Consult conda documentation for further instructions.
 
2. Install BALSAMIC using ``pip`` within the newly created environment: ``pip install -r requirements.txt -e .``

3. Pull container using Singularity: ``singularity pull path_container_file docker://hassanf/balsamic``


If you'd like to install release 5.0.0 the instruction will look like below:

::

    # Checkout a specific tag or branch. If you'd like to install latest changes, use master branch
    git checkout v5.0.0

    # Create a conda env: balsamic_base
    conda env create --file BALSAMIC/conda/balsamic.yaml --name balsamic_base

    # Activate conda environment
    source activate balsamic_base

    # If you don't want to install in editable mode, remove `-e`
    pip install -r requirements.txt -e .

    # Pull container for release_v5.0.0
    singularity pull balsamic_release_v5.0.0 docker://hassanf/balsamic:release_v5.0.0


Automatic Installation
~~~~~~~~~~~~~~~~~~~~~~3

NOTE: The following instructions are for internal use only.  

Use ``install.sh`` script, assuming `${CONDA_ENVS_PATH}` is set to the path for conda environment:

::
  
  ./install.sh -s S -p path_to_conda_env_location -P path_to_container_store_location -c

In above example, the final conda environment will be named: `S_BALSAMIC`

::

    USAGE: ../install.sh [-s _condaprefix -v _balsamic_ver -p _condapath -c]
      1. Conda naming convention: [P,D,S]_[ENVNAME]_%DATE. P: Production, D: Development, S: Stage
      2. Conda environment prefix: Path to conda env. e.g. /home/user/conda_env/
      
      -s _condaprefix     Conda env name prefix. This will be P or D in the help above
      -v _balsamic_ver    Balsamic version tag to install (4.0.0+), or it could be the branch name
      -e _envsuffix       Balsamic conda env suffix. This will be added to the conda env name
      -p _condapath       Conda env path prefix. See point 2 in help above
      -P _containerpath   Container path to store container files. Default set to current directory
      -c                  If set it will use Singularity container for conda instead 
