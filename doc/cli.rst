========
BALSAMIC
========

Bioinformatic Analysis pipeLine for SomAtic MutatIons in Cancer
(**version** = 2.9.6)

.. contents::

BALSAMIC commands
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

    BALSAMIC 2.9.6: Bioinformatic Analysis pipeLine for SomAtic MutatIons in
    Cancer

  Options:
    --version  Show the version and exit.
    --help     Show this message and exit.

  Commands:
    config   create config files required for running the...
    install  Installs required conda environments
    report   Report generator for workflow results
    run      Run BALSAMIC on a provided config file

create config for sample analysis
~~~~~~~~~~~~~

::

    Usage: balsamic config [OPTIONS] COMMAND [ARGS]...

      create config files required for running the pipeline and reporting it

    Options:
      --help  Show this message and exit.

    Commands:
      report  Create a report config file for report generation.
      sample  Create a sample config file from input sample data

::

    Usage: balsamic config sample [OPTIONS]

          Prepares a config file for balsamic run_analysis. For now it is just
          treating json as dictionary and merging them as it is. So this is just
          a placeholder for future.



    Options:
      -a, --analysis-type [paired|single]
                                      Analysis config file for paired (tumor vs
                                      normal) or single (tumor-only) mode.
                                      [default: paired]
      -i, --install-config PATH       Installation config file.  [default: /home/h
                                      assan.foroughi/repo/BALSAMIC/BALSAMIC/config
                                      /install.json]
      -r, --reference-config PATH     Reference config file.  [default: /home/hass
                                      an.foroughi/repo/BALSAMIC/BALSAMIC/config/re
                                      ference.json]
      -p, --panel-bed PATH            Panel bed file for variant calling.
                                      [required]
      -s, --sample-config PATH        Input sample config file.
      -o, --output-config PATH        Output a json config file ready to be
                                      imported for run-analysis  [required]
      -t, --tumor TEXT                Fastq files for tumor sample. Example:
                                      --tumor tumor_1.fastq.gz tumor_2.fastq.gz
      -n, --normal TEXT               Fastq files for normal sample. Example:
                                      --normal normal_1.fastq.gz normal_2.fastq.gz
      --sample-id TEXT                Sample id that is used for reporting, naming
                                      the analysis jobs, and analysis path
      --analysis-dir PATH             Root analysis path to store analysis logs
                                      and results. The final path will be
                                      analysis-dir/sample-id
      --fastq-path PATH               Path for fastq files. All fastq files should
                                      be within same path and that path has to
                                      exist.
      --help                          Show this message and exit.


::

    Usage: balsamic config report [OPTIONS]

      Prepares a config file for balsamic config report to export results as pdf

    Options:
      -s, --sample-config PATH  Sample json config file.  [required]
      -o, --output-config PATH  Path to output config file to write.  [required]
      --help                    Show this message and exit.v

run analysis
~~~~~~~~~~~~

::

    Usage: balsamic run [OPTIONS]

      Runs BALSAMIC workflow on the provided sample's config file

    Options:
      -S, --snake-file PATH      Snakefile required for snakemake to function.
      -s, --sample-config PATH   Sample json config file.  [required]
      -c, --cluster-config PATH  SLURM config json file.  [default: /home/hassan.f
                                 oroughi/repo/BALSAMIC/BALSAMIC/config/cluster.jso
                                 n]
      -l, --log-file PATH        Log file output for BALSAMIC. This is raw log
                                 output from snakemake.
      -r, --run-analysis         By default balsamic run_analysis will run in dry
                                 run mode. Raise thise flag to make the actual
                                 analysis  [default: False]
      -f, --force-all            Force run all analysis. This is same as snakemake
                                 --forceall  [default: False]
      --snakemake-opt TEXT       Pass these options directly to snakemake
      --help                     Show this message and exit.

report generation
~~~~~~~~~~~~

::
    
    Usage: balsamic report [OPTIONS]

    Options:
      -j, --json-report PATH     Input JSON file from workflow output  [required]
      -c, --json-varreport PATH  Input JSON file for variant filters  [required]
      -r, --rulegraph-img PATH   Input rulegraph from workflow output
      --help                     Show this message and exit.

Misc. internal commands
~~~~~~~~~

::

    Usage: balsamic install [OPTIONS]

      Installs conda environments from a conda yaml file.

      By default it doesn't overwrite if the environment by the same name
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



Running variant calling workflow
--------------------------------

In order to run variant calling workflow, first a configuration file

.. code-block:: shell

  balsamic config sample \
    -p path_to_panel_bedfile \
    --sample-id sample_id \
    --normal prefix_to_normal_sample_fastq \
    --tumor prefix_to_tumor_sample_fastq \
    --fastq-path fastq_file_directory \
    --analysis-dir analysis_directory \
    --analysis-type paired_or_single \
    --output-config sample_analysis_config_file_name

The final config file is then set as input for ``run`` subcommand.

.. code-block:: shell

  balsamic run \
    --sample-config sample_analysis_config_file_name -r

After the analysis is finished, the following commands will generate a PDF report

.. code-block:: shell

  balsamic config report --sample-config sample_analysis_config_file_name \
    --output-config sample_report_config_file
    
  balsamic report \
    --json-report sample_report_config_file \
    --json-varreport path_to_BALSAMIC_repo/BALSAMIC/BALSAMIC/config/MSK_impact.json \
 
Config files
------------

BALSAMIC requires two config files: job submission configuration and
analysis configuration. Configurations and their template can be found
within ``config`` directory.
