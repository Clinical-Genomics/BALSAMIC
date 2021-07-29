============
Installation
============

This section describes steps to install BALSAMIC (**version** = 4.1.0)



Software Requirements
~~~~~~~~~~~~~~~~~~~~

- Conda >=version 4.5.0: For detailed software and python requirements please see ``setup.py`` and ``BALSAMIC/conda/balsamic.yaml``
- Singularity >=version 3.0.0: BALSAMIC uses singularity to run vairous parts of the workflow. 
- Python 3.6
- Balsamic is dependent on third-party bioinformatics software ``Sentieon-tools``. Example: for running wgs variant calling using ``TNScope``, and to execute ``UMIworkflow``.

``Note: Set Sentieon envionment variables in your ~/.bashrc file by adding following two lines``

::

    export SENTIEON_INSTALL_DIR=path_to_sentieon_install_dir
    export SENTIEON_LICENSE=IP:Port

Step 1. Installing BALSAMIC
~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Create a conda environment:

::

    conda create -c conda-forge -c defaults --name S_BALSAMIC python==3.7 pip


2. Activate environment: 

::

    conda activate S_BALSAMIC



3. Install BALSAMIC using ``pip`` within the newly created environment: 

::

  pip install --no-cache-dir -U git+https://github.com/Clinical-Genomics/BALSAMIC


Or if you have repository cloned and want it in editable mode:

::

  pip install -e .


Step 2. generate BALSAMIC cache and pull containers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First generate your own COSMIC database key via: https://cancer.sanger.ac.uk/cosmic/help/file_download
The following commands will create and download reference directory at ``~/balsamic_cache`` (change this path if you
want it to be created in another location):

NOTE: This process can take couple of hours

::

  # Note:
  # 1. COSMIC key is in variable $COSMIC_KEY
  # 2. For genome version hg38, set --genome-version to hg38

  balsamic init --outdir ~/balsamic_cache \
    --cosmic-key "${COSMIC_KEY}" \
    --genome-version hg19 \
    --run-analysis
