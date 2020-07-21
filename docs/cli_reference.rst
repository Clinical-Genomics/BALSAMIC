=============
CLI reference
=============

BALSAMIC (**version** = 4.5.0) is using Click. This section covers the CLI reference for BALSAMIC.

.. contents::

Short summary
-------------
Base command, ``balsamic`` has three subcommands: 1) ``install_env`` which is used for installting
conda environemnts 2) ``create_config`` is to create a config file
necessary for running the analysis. 3) ``run_analysis`` is for running
the actual workflow.

Base command
~~~~~~~~~~~~

::

  Usage: balsamic [OPTIONS] COMMAND [ARGS]...

    BALSAMIC 4.5.0: Bioinformatic Analysis pipeLine for SomAtic MutatIons in
    Cancer

  Options:
    --loglevel [DEBUG|INFO|WARNING|ERROR|CRITICAL]
                                    Set the level of log output.  [default:
                                    DEBUG]
    --version                       Show the version and exit.
    --help                          Show this message and exit.

  Commands:
    config   create config files required for running the pipeline.
    plugins  Additional and helper utilities for third party applications
    run      Run BALSAMIC on a provided config file


create config for case analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

  Usage: balsamic config [OPTIONS] COMMAND [ARGS]...

    create config files required for running the pipeline.

  Options:
    --help  Show this message and exit.

  Commands:
    case       Create a sample config file from input sample data
    reference  config workflow for generate reference

::

  Usage: balsamic config case [OPTIONS]

    Prepares a config file for balsamic run_analysis. For now it is just
    treating json as dictionary and merging them as it is. So this is just a
    placeholder for future.

  Options:
    --umi / --no-umi                UMI processing steps for samples with umi
                                    tags  [default: True]
    --umi-trim-length INTEGER       Trim N bases from reads in fastq  [default:
                                    5]
    --quality-trim / --no-quality-trim
                                    Trim low quality reads in fastq  [default:
                                    True]
    --adapter-trim / --no-adapter-trim
                                    Trim adapters from reads in fastq  [default:
                                    False]
    -r, --reference-config PATH     Reference config file.  [required]
    -p, --panel-bed PATH            Panel bed file for variant calling.
    -o, --output-config TEXT        Output a json config filename ready to be
                                    imported for run-analysis
    -t, --tumor TEXT                Fastq files for tumor sample.
                                    Example: if files are
                                    tumor_fqreads_1.fastq.gz
                                    tumor_fqreads_2.fastq.gz,               the
                                    input should be --tumor tumor_fqreads
                                    [required]
    -n, --normal TEXT               Fastq files for normal sample.
                                    Example: if files are
                                    normal_fqreads_1.fastq.gz
                                    normal_fqreads_2.fastq.gz,               the
                                    input should be --normal normal_fqreads
    --case-id TEXT                  Sample id that is used for reporting,
                                    naming the analysis jobs, and analysis path
                                    [required]
    --fastq-prefix TEXT             Prefix to fastq file.               The
                                    string that comes after readprefix
    --analysis-dir PATH             Root analysis path to store
                                    analysis logs and results. The final path
                                    will be analysis-dir/sample-id
    --overwrite-config / --no-overwrite-config
                                    Overwrite output config file
    --create-dir / --no-create-dir  Create analysis directiry.
    --singularity PATH              Download singularity image for BALSAMIC
                                    [required]
    --help                          Show this message and exit.

::

  Usage: balsamic config reference [OPTIONS]

    Configure workflow for reference generation

  Options:
    -o, --outdir TEXT      output directory for ref files eg: reference
                           [required]
    -c, --cosmic-key TEXT  cosmic db authentication key  [required]
    -s, --snakefile PATH   snakefile for reference generation  [default: /home/h
                           assan.foroughi/repos/BALSAMIC/BALSAMIC/workflows/Gene
                           rateRef]
    -d, --dagfile TEXT     DAG file for overview  [default:
                           generate_ref_worflow_graph]
    --singularity PATH     Download singularity image for BALSAMIC  [required]
    --help                 Show this message and exit.

run case analysis and reference creation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

  Usage: balsamic run [OPTIONS] COMMAND [ARGS]...

    Run BALSAMIC on a provided config file

  Options:
    --help  Show this message and exit.

  Commands:
    analysis   Run the analysis on a provided sample config-file
    reference  Run the GenerateRef workflow

::

  Usage: balsamic run analysis [OPTIONS]

    Runs BALSAMIC workflow on the provided sample's config file

  Options:
    -a, --analysis-type [qc|paired|single]
                                    Type of analysis to run from input config
                                    file.              By default it will read
                                    from config file, but it will override
                                    config file               if it is set here.
    -S, --snake-file PATH           Input for a custom snakefile. WARNING: This
                                    is for internal testing,              and
                                    should not be used. Providing a snakefile
                                    supersedes analysis_type option.
    -s, --sample-config PATH        Sample json config file.  [required]
    --run-mode [local|slurm]        Run mode to use. By default SLURM will be
                                    used to run the analysis.              But
                                    local runner also available for local
                                    computing  [default: slurm]
    -c, --cluster-config PATH       SLURM config json file.  [default: /home/has
                                    san.foroughi/repos/BALSAMIC/BALSAMIC/config/
                                    cluster.json]
    -l, --log-file PATH             Log file output for BALSAMIC.
                                    This is raw log output from snakemake.
    -r, --run-analysis              By default balsamic run_analysis will run in
                                    dry run mode.               Raise thise flag
                                    to make the actual analysis  [default:
                                    False]
    --qos [low|normal|high]         QOS for sbatch jobs. Passed to /home/hassan.
                                    foroughi/repos/BALSAMIC/BALSAMIC/commands/ru
                                    n/sbatch.py  [default: low]
    -f, --force-all                 Force run all analysis. This is same as
                                    snakemake --forceall  [default: False]
    --snakemake-opt TEXT            Pass these options directly to snakemake
    --slurm-account TEXT            SLURM account to run jobs
    --slurm-mail-user TEXT          SLURM mail user to send out email.
    --slurm-mail-type [NONE|BEGIN|END|FAIL|REQUEUE|ALL|TIME_LIMIT]
                                    SLURM mail type to send out email.
                                    This will be applied to all jobs and
                                    override snakemake settings.
    --help                          Show this message and exit.

::

  Usage: balsamic run reference [OPTIONS]

    Run generate reference workflow

  Options:
    -s, --snakefile TEXT      snakefile for reference generation
    -c, --configfile TEXT     Config file to run the workflow  [required]
    --run-mode [slurm|local]  Run mode to use.(LOCAL, SLURM for HPC)
    --cluster-config PATH     SLURM config json file.  [default: /home/hassan.fo
                              roughi/repos/BALSAMIC/BALSAMIC/config/cluster.json
                              ]
    -l, --log-file PATH       Log file output for BALSAMIC. This is raw log
                              output from snakemake.
    -r, --run-analysis        By default balsamic run_analysis will run in dry
                              run mode.               Raise thise flag to make
                              the actual analysis  [default: False]
    --qos [low|normal|high]   QOS for sbatch jobs. Passed to /home/hassan.foroug
                              hi/repos/BALSAMIC/BALSAMIC/commands/run/sbatch.py
                              [default: low]
    -f, --force-all           Force run all analysis. This is same as snakemake
                              --forceall  [default: False]
    --snakemake-opt TEXT      Pass these options directly to snakemake
    --help                    Show this message and exit.
