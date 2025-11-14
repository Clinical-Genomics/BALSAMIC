============
Installation
============

This section describes steps to install BALSAMIC (**version** = 18.0.0)



Software Requirements
~~~~~~~~~~~~~~~~~~~~~

- Conda >=version 4.5.0: For detailed software and python requirements please see ``setup.py`` and ``BALSAMIC/conda/balsamic.yaml``
- Singularity >=version 3.0.0: BALSAMIC uses singularity to run various parts of the workflow.
- Python 3.11
- BALSAMIC is dependent on third-party bioinformatics software ``Sentieon-tools`` for all workflows.
- The BALSAMIC wrapper is hard-coded to the SLURM workload manager and requires scontrol

``Note: To run Balsamic you need to supply the --sentieon-install-dir and --sentieon-license arguments during the config``


Step 1. Installing BALSAMIC
~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Create a conda environment:

::

    conda create --name S_balsamic -c conda-forge python=3.11 pip "cython<3" pygraphviz wkhtmltopdf snakemake-executor-plugin-slurm=1.9.2



2. Activate environment:

::

    conda activate S_BALSAMIC



3. Install BALSAMIC using ``pip`` within the newly created environment:

::

  pip install --no-build-isolation --no-cache-dir -U git+https://github.com/Clinical-Genomics/BALSAMIC


Or if you have repository cloned and want it in editable mode:

::

  pip install -e .


Step 2. generate BALSAMIC cache and pull containers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First generate your own COSMIC database key via: https://cancer.sanger.ac.uk/cosmic/register
The following commands will create and download reference directory at ``~/balsamic_cache`` (change this path if you
want it to be created in another location):

NOTE: This process can take couple of hours

::

  # Note:
  # 1. COSMIC key is in variable $COSMIC_KEY
  # 2. For genome version hg38, set --genome-version to hg38
  # 3. For using develop container version, set --cache-version to develop
  # 4. For submitting jobs to slurm cluster, use option --account

  balsamic init --outdir ~/balsamic_cache \
    --cosmic-key "${COSMIC_KEY}" \
    --genome-version hg19 \
    --run-analysis \
    --account development

  # Generate cache locally instead of slurm job submission
  balsamic init --outdir ~/balsamic_cache \
    --cosmic-key "${COSMIC_KEY}" \
    --genome-version hg19 \
    --run-analysis \
    --run-mode local \
    --snakemake-opt "--cores 16"

