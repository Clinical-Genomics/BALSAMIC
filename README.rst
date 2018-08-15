========
BALSAMIC
========

Bioinformatic Analysis pipeLine for SomAtic MutatIons in Cancer
(**version** = 1.13.0)

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

Usage
-----

After a successfull installation, BALSAMIC is within
``D_BALSAMIC_${prefix}`` conda environment. Base command, ``balsamic``
has three subcommands: 1) ``install_env`` which is used for installting
conda environemnts 2) ``create_config`` is to create a config file
necessary for running the analysis. 3) ``run_analysis`` is for running
the actual workflow.

Base command
~~~~~~~~~~~~

::

   Usage: balsamic [OPTIONS] COMMAND [ARGS]...

     BALSAMIC 1.0.3-rc2: Bioinformatic Analysis pipeLine for SomAtic MutatIons
     in Cancer

   Options:
     --version  Show the version and exit.
     --help     Show this message and exit.

   Commands:
     create_config  Create a sample config file from input sample data
     install_env    Installs required conda environments
     run_analysis   Run BALSAMIC on a provided config file

install_env
~~~~~~~~~~~

::

   Usage: balsamic install_env [OPTIONS]

     Installs conda environments from a conda yaml file.

     By default it doesn''t overwrite if the environment by the same name
     exists. If _overwrite_ flag is provided, it tries to remove the enviroment
     first, and then install it in the path provided.

   Options:
     -i, --input-conda-yaml PATH     Input conda yaml file.  [required]
     -s, --env-name-suffix TEXT      Mandatory alphanumeric suffix for
                                     environment name.  [required]
     -o, --overwrite-env             Overwite conda enviroment if it exists.
                                     Default = no. WARNING: The environment with
                                     matching name will be deleted  [default:
                                     False]
     -d, --env-dir-prefix TEXT       Conda enviroment directory. It will be
                                     ignored if its provided within yaml file.
                                     Format: /path/env/envname.
     -p, --packages-output-yaml PATH
                                     Output a yaml file containing packages
                                     installed in each input yaml file.
                                     [required]
     --help                          Show this message and exit.

create_config
~~~~~~~~~~~~~

::

   Usage: balsamic create_config [OPTIONS]

         Prepares a config file for balsamic run_analysis. For now it is just
         treating json as dictionary and merging them as it is. So this is just
         a placeholder for future.



   Options:
     -a, --analysis-config PATH   Analysis config file.  [required]
     -i, --install-config PATH    Installation config file.  [required]
     -r, --reference-config PATH  Reference config file.  [required]
     -s, --sample-config PATH     Input sample config file.  [required]
     -o, --output-config PATH     Output a json config file ready to be imported
                                  for run-analysis  [required]
     --help                       Show this message and exit.

run_analysis
~~~~~~~~~~~~

::

    Usage: balsamic run_analysis [OPTIONS]
    
        Runs BALSAMIC workflow on the provided sample''s config file
    
    Options:
    -S, --snake-file PATH      Snakefile required for snakemake to function.
                            [required]
    -s, --sample-config PATH   Sample json config file.  [required]
    -c, --cluster-config PATH  SLURM config json file.  [required]
    -l, --log-file PATH        Log file output for BALSAMIC. This is raw log
                            output from snakemake.
    -r, --run-analysis         By default balsamic run_analysis will run in dry
                            run mode. Raise thise flag to make the actual
                            analysis  [default: False]
    -f, --force-all            Force run all analysis. This is same as snakemake
                            --forceall  [default: False]
    --snakemake-opt TEXT       Pass these options directly to snakemake
    --help                     Show this message and exit.
    

Running variant calling workflow
--------------------------------

In order to run variant calling workflow, first a configuration file
must be created. It requires a sample.json file, and panel data file. A template for sample.json can be found within
`BALSAMIC/config/sample.json`. Otherwise, the following parameters must be provided. A example with or without
sample.json file is as below

.. code-block:: shell

   balsamic create_config \
     --panel-bed path_to_panel_bed_file \
     --sample-config BALSAMIC/config/sample.json \
     --output-config BALSAMIC/config/sample_analysis.json

.. code-block:: shell
    
    balsamic create_config \
      --normal base_name_to_normal_sample \
      --tumor base_name_to_tumor_sample \
      --sample-id sample_name \
      --analysis_type paired \
      --analysis-dir path_to_store_analysis_dir \
      --fastq-path path_to_fastq_files_for_tumor_and_normal \
      --output-config path_and_filename_output_config_file \
      --panel-bed path_to_panel_bed_file

The final config file is then set as input for ``run_analysis``
subcommand.

.. code-block:: shell

   balsamic run_analysis \
     --sample-config BALSAMIC/config/sample_analysis.json

Config files
------------

BALSAMIC requires two config files: job submission configuration and
analysis configuration. Configurations and their template can be found
within ``config`` directory. The only config file that user needs to
provide is the ``sample.json`` by updating the necessary entries.
