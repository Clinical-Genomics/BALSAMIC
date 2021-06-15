[7.X.X]
-------

Added:
^^^^^^

* samtools flagstats and stats to workflow and MultiQC
* Delly v0.8.7 somatic SV caller
* Delly containter
* bcftools v1.12 to delly container
* tabix v0.2.6 to delly container
* PASSed SV calls from Manta to clinical delivery
* An extra filter to VarDict tumor-normal to remove variants with STATUS=Germline, all other will still be around
* Added ``vcf2cytosure``` to annotate container
* ``git`` to the container definition

Changed:
^^^^^^^^

* Upgrade to latest sentieon version 202010.02
* New name ``MarkDuplicates`` to ``picard_markduplicates`` in ``bwa_mem`` rule and ``cluster.json``
* New name rule ``GATK_contest`` to ``gatk_contest`` 
* Avoid running pytest github actions workflow on ``docs/**`` and ``CHANGELOG.rst`` changes

Fixed:
^^^^^^

* post-processing of the umi consensus in handling BI tags
* vcf-filtered-clinical tag files will have all variants including PASS
* Refactor snakemake align rules according to snakemake etiquette 
* Refactor snakemake fastqc vep contest and mosdepth rules according to snakemake etiquette
* Refactor snakemake manta rule according to snakemake etiquette

Removed:
^^^^^^^^

* Cleaned up unused container definitions and conda environment files

[7.2.2]
-------

Fixed:
^^^^^^

* An error with Sentieon for better management of memory fixes #621

[7.2.1]
-------

Changed:
^^^^^^^^

* Rename Github actions to reflect their content

[7.2.0]
-------

Added:
^^^^^^

* Changelog reminder workflow to Github
* Snakemake workflow for created PON reference
* Balsamic cli config command(pon) for creating json for PON analysis
* tumor lod option for passing tnscope-umi final variants
* Git guide to make balsamic release in FAQ docs

Changed:
^^^^^^^^

* Expanded multiqc result search dir to whole analysis dir
* Simple test for docker container

Fixed:
^^^^^^

* Correctly version bump for Dockerfile

Removed:
^^^^^^^^

* Removed unused Dockerfile releases
* Removed redundant genome version from ``reference.json``

[7.1.10]
-------

Fixed:
^^^^^^

* Bug in ``ngs_filter`` rule set for tumor-only WGS
* Missing delivery of tumor only WGS filter

[7.1.9]
-------


Changed:
^^^^^^^^

* only pass variants are not part of delivery anymore
* delivery tag file ids are properly matched with sample_name
* tabix updated to 0.2.6
* fastp updated to 0.20.1
* samtools updated to 1.12
* bedtools updated to 2.30.0

Removed:
^^^^^^^^

* sentieon-dedup rule from delivery
* Removed all pre filter pass from delivery


[7.1.8]
-------

Fixed:
^^^^^^

* Target coverage (Picard HsMetrics) for UMI files is now correctly calculated.

Changed:
^^^^^^^^


* TNscope calculated AF values are fetched and written to AFtable.txt.

[7.1.7]
-------

Added:
^^^^^^

* ngs_filter_tnscope is also part of deliveries now

Changed:
^^^^^^^^

* rankscore is now a research tag instead of clinical
* Some typo and fixes in the coverage and constant metrics
* Delivery process is more verbose

Fixed:
^^^^^^

* CNVKit output is now properly imported in the deliveries and workflow

[7.1.6]
-------

Fixed:
^^^^^^

* CSS style for qc coverage report is changed to landscape

[7.1.5]
-------

Changed:
^^^^^^^^

* update download url for 1000genome WGS sites from ftp to http

[7.1.4]
-------

Changed:
^^^^^^^^

* bump picard to version 2.25.0

[7.1.3]
-------

Fixed:
^^^^^

* ``assets`` path is now added to bind path

[7.1.2]
-------

Fixed:
^^^^^

* umi_workflow config json is set as true for panel and wgs as false.
* Rename umiconsensus bam file headers from {samplenames} to TUMOR/NORMAL. 
* Documentation autobuild on RTFD


[7.1.1]
-------

Fixed:
^^^^^

* Moved all requirements to setup.py, and added all package_data there. Clean up unused files.

[7.1.0]
-------

Removed
^^^^^^^

* ``tnsnv`` removed from WGS analysis, both tumor-only and tumor-normal
* GATK-BaseRecalibrator is removed from all workflows

Fixed
^^^^^

* Fixed issue 577 with missing ``tumor.merged.bam`` and ``normal.merged.bam`` 
* Issue 448 with lingering tmp_dir. It is not deleted after analysis is properly finished.

Changed
^^^^^^^

* All variant calling rules use proper ``tumor.merged.bam`` or ``normal.merged.bam`` as inputs

[7.0.2]
-------

Added
^^^^^

* Updated docs with FAQ for UMI workflow

Fixed
^^^^^

* fix job scheduling bug for benchmarking
* rankscore's output is now a proper vcf.gz file
* Manta rules now properly make a sample_name file


[7.0.1]
-------

Added
^^^^^

* github action workflow to autobuild release containers


[7.0.0]
-------

Added
^^^^^

* ``balsamic init`` to download reference and related containers done in PRs #464 #538
* ``balsamic config case`` now only take a cache path instead of container and reference #538
* UMI workflow added to main workflow in series of PRs #469 #477 #483 #498 #503 #514 #517
* DRAGEN for WGS applications in PR #488
* A framework for QC check PR #401
* ``--quiet``` option for ``run analysis`` PR #491
* Benchmark SLURM jobs after the analysis is finished PR #534
* One container per conda environment (i.e. decouple containers) PR #511 #525 #522
* ``--disable-variant-caller`` command for ``report deliver`` PR #439
* Added genmod and rankscore in series of two PRs #531 and #533
* Variant filtering to Tumor-Normal in PR #534
* Split SNV/InDels and SVs from TNScope variant caller PR #540
* WGS Tumor only variant filters added in PR #548

Changed
^^^^^^^

* Update Manta to 1.6.0 PR #470
* Update FastQC to 0.11.9 PR #532
* Update BCFTools to 1.11 PR #537
* Update Samtools to 1.11 PR #537
* Increase resources and runtime for various workflows in PRs #482 
* Python package dependenicies versions fixed in PR #480
* QoL changes to workflow in series of PR #471
* Series of documentation updates in PRs #489 #553
* QoL changes to scheduler script PR #491
* QoL changes to how temporary directories are handlded PR #516
* TNScope model apply rule merged with TNScope variant calling for tumor-normal in WGS #540
* Decoupled ``fastp`` rule into two rules to make it possible to use it for UMI runs #570


Fixed
^^^^^

* A bug in Manta variant calling rules that didn't name samples properly to TUMOR/NORMAL in the VCF file #572


[6.1.2]
-------

Changed
^^^^^^^
* Changed hk delivery tag for coverage-qc-report


[6.1.1]
-------

Fixed
^^^^^

* No UMI trimming for WGS applications #486
* Fixed a bug where BALSAMIC was checking for sacct/jobid file in local mode PR #497
* ``readlink`` command in ``vep_germline``, ``vep_somatic``, ``split_bed``, and ``GATK_popVCF`` #533
* Fix various bugs for memory handling of Picardtools and its executable in PR #534
* Fixed various issues with ``gsutils`` in PR #550

Removed
^^^^^^^

* ``gatk-register`` command removed from installing GATK PR #496

[6.1.1]
-------

* Fixed a bug with missing QC templates after ``pip install``


[6.1.0]
-------

Added
^^^^^
* CLI option to expand report generation for TGA and WES runs. Please see ``balsamic report deliver --help``
* BALSAMIC now generates a custom HTML report for TGA and WES cases.


[6.0.4]
-------

Changed
^^^^^^^

* Reduces MQ cutoff from 50 to 40 to only remove obvious artifacts PR #535
* Reduces AF cutoff from 0.02 to 0.01 PR #535

[6.0.3]
-------

Added
^^^^^

* ``config case`` subcommand now has ``--tumor-sample-name`` and ``--normal-sample-name``

Fixed
^^^^^

* Manta resource allocation is now properly set PR #523
* VarDict resource allocation in cluster.json increased (both core and time allocation) PR #523
* minimum memory request for GATK mutect2 and haplotypecaller is removed and max memory increased PR #523

[6.0.2]
-------

Added
^^^^^

* Document for Snakemake rule grammar PR #489


Fixed
^^^^^

* removed ``gatk3-register`` command from Dockerfile(s) PR #508


[6.0.1]
-------

Added
^^^^^
* A secondary path for latest jobids submitted to cluster (slurm and qsub) PR #465

[6.0.0]
-------

Added
^^^^^
* UMI workflow using Sentieon tools. Analysis run available via `balsamic run analysis --help` command. PR #359
* VCFutils to create VCF from flat text file. This is for internal purpose to generate validation VCF. PR #349
* Download option for hg38 (not validated) PR #407
* Option to disable variant callers for WES runs. PR #417

Fixed
^^^^^
* Missing cyvcf2 dependency, and changed conda environment for base environment PR #413
* Missing numpy dependency PR #426

Changed
^^^^^^^
* COSMIC db for hg19 updated to v90 PR #407
* Fastp trimming is now a two-pass trimming and adapter trimming is always enabled. This might affect coverage slightly PR #422
* All containers start with a clean environment #425
* All Sentieon environment variables are now added to config when workflow executes #425
* Branching model will be changed to gitflow

[5.1.0]
-------

Fixed
^^^^^
* Vardict-java version fixed. This is due to bad dependency and releases available on conda. Anaconda is not yet update with vardict 1.8, but vardict-java 1.8 is there. This causes various random breaks with Vardict's TSV output. #403

Changed
^^^^^^^
* Refactored Docker files a bit, preparation for decoupling #403

Removed
^^^^^^^
* In preparation for GATK4, IndelRealigner is removed #404


[5.0.1]
-------

Added
^^^^^
* Temp directory for various rules and workflow wide temp directory #396

Changed
^^^^^^^
* Refactored tags for housekeeper delivery to make them unique #395
* Increased core requirements for mutect2 #396
* GATK3.8 related utils run via jar file instead of gatk3 #396


[5.0.0]
-------

Added
^^^^^
* Config.json and DAG draph included in Housekeeper report #372
* New output names added to cnvkit_single and cnvkit_paired #372
* New output names added to vep.rule #372
* Delivery option to CLI and what to delivery with delivery params in rules that are needed to be delivered #376
* Reference data model with validation #371
* Added container path to install script #388

Changed
^^^^^^^
* Delivery file format simplified #376
* VEP rules have "all" and "pass" as output #376
* Downloaded reference structure changed #371
* genome/refseq.flat renamed to genome/refGene.flat #371
* reverted CNVKit to version 0.9.4 #390

Fixed
^^^^^
* Missing pygments to requirements.txt to fix travis CI #364
* Wildcard resolve for deliveries of vep_germline #374
* Missing index file from deliverables #383
* Ambiguous deliveries in vep_somatic and ngs_filters #387
* Updated documentation to match with installation #391

Removed
^^^^^^^
* Temp files removed from list of outputs in vep.rule #372
* samtools.rule and merged it with bwa_mem #375


[4.5.0]
-------

Added
^^^^^
* Models to build config case JSON. The models and descriptions of their contents can now be found
  in BALSAMIC/utils/models.py
* Added analysis_type to `report deliver` command
* Added report and delivery capability to Alignment workflow
* run_validate.sh now has -d to handle path to analysis_dir (for internal use only) #361

Changed
^^^^^^^

* Fastq files are no longer being copied as part of creation of the case config file.
  A symlink is now created at the destination path instead
* Config structure is no longer contained in a collestion of JSON files.
  The config models are now built using Pydantic and are contained in BALSAMIC/utils/models.py

Removed
^^^^^^^

* Removed command line option "--fastq-prefix" from config case command
* Removed command line option "--config-path" from config case command.
  The config is now always saved with default name "case_id.json"
* Removed command line option "--overwrite-config" from config-case command
  The command is now always executed with "--overwrite-config True" behavior

Refactored
^^^^^^^^^^

* Refactored BALSAMIC/commands/config/case.py:
  Utility functions are moved to BALSAMIC/utils/cli.py
  Models for config fields can be found at BALSAMIC/utils/models.py
  Context aborts and logging now contained in pilot function
  Tests created to support new architecture
* Reduce analysis directory's storage

Fixed
^^^^^
* Report generation warnings supressed by adding workdirectory
* Missing tag name for germline annotated calls #356
* Bind path is not added as None if analysis type is wgs #357
* Changes vardict to vardict-java #361


[4.4.0]
-------

Added
^^^^^

* pydantic to validate various models namely variant caller filters

Changed
^^^^^^^

* Variant caller filters moved into pydantic
* Install script and setup.py
* refactored install script with more log output and added a conda env suffix option
* refactored docker container and decoupled various parts of the workflow


[4.3.0]
-------


Added
^^^^^

* Added cram files for targeted sequencing runs fixes #286
* Added `mosdepth` to calculate coverage for whole exome and targeted sequencing
* Filter models added for tumor-only mode
* Enabling adapter trim enables pe adapter trim option for fastp
* Annotate germline variant calls
* Baitset name to picard hsmetrics

Deprecated
^^^^^^^^^^

* Sambamba coverage and rules will be deprecated

Fixed
^^^^^

* Fixed latest tag in install script
* Fixed lack of naming final annotated VCF TUMOR/NORMAL


Changed
^^^^^^^

* Increased run time for various slurm jobs fixes #314
* Enabled SV calls for VarDict tumor-only
* Updated `ensembl-vep` to v100.2

[4.2.4]
-------


Fixed
^^^^^

* Fixed sort issue with bedfiles after 100 slop


[4.2.3]
-------

Added
^^^^^


* Added Docker container definition for release and bumpversion

Changed
^^^^^^^


* Quality of life change to rtfd docs

Fixed
^^^^^


* Fix Docker container with faulty git checkout

[4.2.2]
-------

Added
^^^^^


* Add "SENTIEON_TMPDIR" to wgs workflow

[4.2.1]
-------

Changed
^^^^^^^


* Add docker container pull for correct version of install script

[4.2.0]
-------

Added
^^^^^


* CNV output as VCF
* Vep output for PASSed variants
* Report command with status and delivery subcommands

Changed
^^^^^^^


* Bed files are slopped 100bp for variant calling fix #262
* Disable vcfmerge
* Picard markduplicate output moved from log to output
* Vep upgraded to 99.1
* Removed SVs from vardict
* Refactored delivery plugins to produce a file with list of output files from workflow
* Updated snakemake to 5.13

Fixed
^^^^^


* Fixed a bug where threads were not sent properly to rules

Removed
^^^^^^^


* Removed coverage annotation from mutect2
* Removed source deactivate from rules to suppress conda warning
* Removed ``plugins delivery`` subcommand
* Removed annotation for germline caller results

[4.1.0]
-------

Added
^^^^^


* VEP now also produces a tab delimited file
* CNVkit rules output genemetrics and gene break file
* Added reference genome to be able to calculate AT/CG dropouts by Picard
* coverage plot plugin part of issue #75
* callable regions for CNV calling of tumor-only

Changed
^^^^^^^


* Increased time for indel realigner and base recalib rules
* decoupled vep stat from vep main rule
* changed qsub command to match UGE
* scout plugin updated

Fixed
^^^^^


* WGS qc rules - updated with correct options
  (picard - CollectMultipleMetrics, sentieon - CoverageMetrics)
* Log warning if WES workflow cannot find SENTIEON* env variables
* Fixes issue with cnvkit and WGS samples #268
* Fix #267 coverage issue with long deletions in vardict

[4.0.1] - 2019-11-08
--------------------

Added
^^^^^


* dependencies for workflow report
* sentieon variant callers germline and somatic for wes cases

Changed
^^^^^^^


* housekeeper file path changed from basename to absolute
* scout template for sample location changed from delivery_report to scout
* rule names added to benchmark files

[4.0.0] - 2019-11-04
--------------------

SGE qsub support release

Added
^^^^^


* ``install.sh`` now also downloads latest container
* Docker image for balsamic as part of ci
* Support for qsub alongside with slurm on ``run analysis --profile``

Changed
^^^^^^^


* Documentation updated
* Test fastq data and test panel bed file with real but dummy data

[3.3.1] - 2019-10-28
--------------------

Fixed
^^^^^


* Various links for reference genome is updated with working URL
* Config reference command now print correct output file

[3.3.0] - 2019-10-24
--------------------

somatic vcfmerge release

Added
^^^^^


* QC metrics for WGS workflow
* refGene.txt download to reference.json and reference workflow
* A new conda environment within container
* A new base container built via Docker (centos7:miniconda3_4_6_14)
* VCFmerge package as VCF merge rule (https://github.com/hassanfa/VCFmerge)
* A container for develop branch
* Benchmark rules to variant callers

Changed
^^^^^^^


* SLURM resource allocation for various variancalling rules optimized
* mergetype rule updated and only accepts one single tumor instead of multiple

[3.2.3] - 2019-10-24
--------------------

Fixed
^^^^^


* Removed unused output files from cnvkit which caused to fail on targetted analysis

[3.2.2] - 2019-10-23
--------------------

Fixed
^^^^^


* Removed target file from cnvkit batch

[3.2.1] - 2019-10-23
--------------------

Fixed
^^^^^


* CNVkit single missing reference file added

[3.2.0] - 2019-10-11
--------------------

Adds:
^^^^^


* CNVkit to WGS workflow
* get_thread for runs

Changed:
^^^^^^^^


* Optimized resources for SLURM jobs

Removed:
^^^^^^^^


* Removed hsmetrics for non-mark duplicate bam files

[3.1.4] - 2019-10-08
--------------------

Fixed
^^^^^


* Fixes a bug where missing capture kit bed file error for WGS cases

[3.1.3] - 2019-10-07
--------------------

Fixed
^^^^^


* benchmark path bug issue #221

[3.1.2] - 2019-10-07
--------------------

Fixed
^^^^^


* libreadline.so.6 symlinking and proper centos version for container

[3.1.1] - 2019-10-03
--------------------

Fixed
^^^^^


* Proper tag retrieval for release
  ### Changed
* BALSAMIC container change to latest and version added to help line

[3.1.0] - 2019-10-03
--------------------

TL;DR:


* QoL changes to WGS workflow
* Simplified installation by moving all tools to a container

Added
^^^^^


* Benchmarking using psutil
* ML variant calling for WGS
* ``--singularity`` option to ``config case`` and ``config reference``

Fixed
^^^^^


* Fixed a bug with boolean values in analysis.json

Changed
^^^^^^^


* ``install.sh`` simplified and will be depricated
* Singularity container updated
* Common somatic and germline variant callers are put in single file
* Variant calling workflow and analysis config files merged together

Removed
^^^^^^^


* ``balsamic install`` is removed
* Conda environments for py36 and py27 are removed

[3.0.1] - 2019-09-11
--------------------

Fixed
^^^^^


* Permissions on ``analysis/qc`` dir are 777 now

[3.0.0] - 2019-09-05
--------------------

This is major release.
TL;DR:


* Major changes to CLI. See documentation for updates.
* New additions to reference generation and reference config file generation and complete overhaul
* Major changes to reposityory structure, conda environments.

Added
^^^^^


* Creating and downloading reference files: ``balsamic config reference`` and ``balsamic run reference``
* Container definitions for install and running BALSAMIC
* Bunch of tests, setup coveralls and travis.
* Added Mutliqc, fastp to rule utilities
* Create Housekeeper and Scout files after analysis completes
* Added Sentieon tumor-normal and tumor only workflows
* Added trimming option while creating workflow
* Added multiple tumor sample QC analysis
* Added pindle for indel variant calling
* Added Analysis finish file in the analysis directory

Fixed
^^^^^


* Multiple fixes to snakemake rules

Changed
^^^^^^^


* Running analysis through: ``balsamic run analysis``
* Cluster account and email info added to ``balsamic run analysis``
* ``umi`` workflow through ``--umi`` tag. [workflow still in evaluation]
* ``sample-id`` replaced by ``case-id``
* Plan to remove FastQC as well

Removed
^^^^^^^


* ``balsamic config report`` and ``balsamic report``
* ``sample.config`` and ``reference.json`` from config directory
* Removed cutadapt from workflows

[2.9.8] - 2019-01-01
--------------------

Fixed
^^^^^


* picard hsmetrics now has 50000 cov max
* cnvkit single wildcard resolve bug fixed

[2.9.7] - 2019-02-28
--------------------

Fixed
^^^^^


* Various fixes to umi_single mode
* analysis_finish file does not block reruns anymore
* Added missing single_umi to analysis workflow cli

Changed
^^^^^^^


* vardict in single mode has lower AF threshold filter (0.005 -> 0.001)

[2.9.6] - 2019-02-25
--------------------

Fixed
^^^^^


* Reference to issue #141, fix for 3 other workflows
* CNVkit rule update for refflat file

[2.9.5] - 2019-02-25
--------------------

Added
^^^^^


* An analysis finish file is generated with date and time inside (%Y-%M-%d T%T %:z)

[2.9.4] - 2019-02-13
--------------------

Fixed
^^^^^


* picard version update to 2.18.11 github.com/hassanfa/picard

[2.9.3] - 2019-02-12
--------------------

Fixed
^^^^^


* Mutect single mode table generation fix
* Vardict single mode MVL annotation fix

[2.9.2] - 2019-02-04
--------------------

Added
^^^^^


* CNVkit single sample mode now in workflow
* MVL list from cheng et al. 2015 moved to assets

[2.9.1] - 2019-01-22
--------------------

Added
^^^^^


* Simple table for somatic variant callers for single sample mode added

Fixed
^^^^^


* Fixes an issue with conda that unset variables threw an error issue #141

[2.9.0] - 2019-01-04
--------------------

Changed
^^^^^^^


* Readme structure and example
* Mutect2's single sample output is similar to paired now
* cli path structure update

Added
^^^^^


* test data and sample inputs
* A dag PDF will be generated when config is made
* umi specific variant calling

[2.8.1] - 2018-11-28
--------------------

Fixed
^^^^^


* VEP's perl module errors
* CoverageRep.R now properly takes protein_coding transcatipts only

[2.8.0] - 2018-11-23
--------------------

UMI single sample align and QC

Added
^^^^^


* Added rules and workflows for UMI analysis: QC and alignment

[2.7.4] - 2018-11-23
--------------------

Germline single sample

Added
^^^^^


* Germline single sample addition
  ### Changed
* Minor fixes to some rules to make them compatible with tumor mode

[2.7.3] - 2018-11-20
--------------------

Fixed
^^^^^


* Various bugs with DAG to keep popvcf and splitbed depending on merge bam file
* install script script fixed and help added

[2.7.2] - 2018-11-15
--------------------

Changed
^^^^^^^


* Vardict, Strelka, and Manta separated from GATK best practice pipeline

[2.7.1] - 2018-11-13
--------------------

Fixed
^^^^^


* minro bugs with strelka_germline and freebayes merge
  ### Changed
* removed ERC from haplotypecaller

[2.7.0] - 2018-11-08
--------------------

Germline patch

Added
^^^^^


* Germline caller tested and added to the paired analysis workflow: Freebayes, HaplotypeCaller, Strelka, Manta

Changed
^^^^^^^


* Analysis config files updated
* Output directory structure changed
* vep rule is now a single rule
* Bunch of rule names updated and shortened, specifically in Picard and GATK
* Variant caller rules are all updated and changed
* output vcf file names are now more sensible: {SNV,SV}.{somatic,germline}.sampleId.variantCaller.vcf.gz
* Job limit increased to 300

Removed
^^^^^^^


* removed bcftools.rule for var id annotation

Changed
^^^^^^^

Fixed
^^^^^

[2.6.3] - 2018-11-01
--------------------

Changed
^^^^^^^


* Ugly and godforsaken ``runSbatch.py`` is now dumping sacct files with job IDs. Yikes!

[2.6.2] - 2018-10-31
--------------------

Fixed
^^^^^


* added ``--fastq-prefix`` option for ``config sample`` to set fastq prefix name. Linking is not changed.

[2.6.1] - 2018-10-29
--------------------

Fixed
^^^^^


* patched a bug for copying results for strelka and manta which was introduced in ``2.5.0``

[2.5.0] - 2018-10-22
--------------------

Changed
^^^^^^^


* ``variant_panel`` changed to ``capture_kit``
* sample config file takes balsamic version
* bioinfo tool config moved bioinfotool to cli_utils from ``config report``

Added
^^^^^


* bioinfo tool versions is now added to analysis config file

[2.4.0] - 2018-10-22
--------------------

Changed
^^^^^^^


* ``balsamic run`` has 3 stop points: paired variant calling, single mode variant calling, and QC/Alignment mode.
* ``balsamic run [OPTIONS] -S ...`` is depricated, but it supersedes ``analysis_type`` mode if provided.

[2.3.3] - 2018-10-22
--------------------

Added
^^^^^


* CSV output for variants in each variant caller based on variant filters
* DAG image of workflow
  ### Changed
* Input for variant filter has a default value
* ``delivery_report`` is no created during config generation
* Variant reporter R script cmd updated in ``balsamic report``

[2.3.2] - 2018-10-19
--------------------

Changed
^^^^^^^


* Fastq files are now always linked to ``fastq`` directory within the analysis directory

Added
^^^^^


* ``balsamic config sample`` now accepts individual files and paths. See README for usage.

[2.3.1] - 2018-09-25
--------------------

Added
^^^^^


* CollectHSmetric now run twice for before and after markduplicate

[2.3.0] - 2018-09-25
--------------------

Changed
^^^^^^^


* Sample config file now includes a list of chromosomes in the panel bed file

Fixed
^^^^^


* Non-matching chrom won't break the splitbed rule anymore
* collectqc rules now properly parse tab delimited metric files

[2.2.0] - 2018-09-11
--------------------

Added
^^^^^


* Coverage plot to report
* target coverage file to report json
* post-cutadapt fastqc to collectqc
* A header to report pdf
* list of bioinfo tools used in the analysis added to report
  ### Changed
* VariantRep.R now accepts multiple inputs for each parameter (see help)
* AF values for MSKIMPACT config
  ### Fixed
* Output figure for coverageplot is now fully square :-)

[2.1.0] - 2018-09-11
--------------------

Added
^^^^^


* normalized coverage plot script
* fastq file IO check for config creation
* added qos option to ``balsamic run``
  ### Fixed
* Sambamba depth coverage parameters
* bug with picard markduplicate flag

[2.0.2] - 2018-09-11
--------------------

Added
^^^^^


* Added qos option for setting qos to run jobs with a default value of low

[2.0.1] - 2018-09-10
--------------------

Fixed
^^^^^


* Fixed package dependencies with vep and installation

[2.0.0] - 2018-09-05
--------------------

Variant reporter patch and cli update

Added
^^^^^


* Added ``balsamic config sample`` and ``balsamic config report`` to generate run analysis and reporting config
* Added ``VariantRep.R`` script to information from merged variant table: variant summry, TMB, and much more
* Added a workflow for single sample mode alignment and QC only
* Added QC skimming script to qccollect to generate nicely formatted information from picard
  ### Changed
* Change to CLI for running and creating config
* Major overhaul to coverage report script. It's now simpler and more readable!
  ### Fixed
* Fixed sambamba depth to include mapping quality
* Markduplicate now is now by default on marking mode, and will NOT remove duplicates
* Minor formatting and script beautification happened

[1.13.1] - 2018-08-17
---------------------

Fixed
^^^^^


* fixed a typo in MSKMVL config
* fixed a bug in strelka_simple for correct column orders

[1.13.0] - 2018-08-10
---------------------

Added
^^^^^


* rule for all three variant callers for paired analysis now generate a simple VCF file
* rule for all three variant callers for paired analysis to convert VCF into table format
* MVL config file and MVL annotation to VCF calls for SNV/INDEL callers
* CALLER annotation added to SNV/INDEL callers
* exome specific option for strelka paired
* create_config subcommand is now more granular, it accepts all enteries from sample.json as commandline arguments
* Added tabQuery to the assets as a tool to query the tabulated output of summarized VCF
* Added MQ annotation field to Mutect2 output see #67
  ### Changed
* Leaner VCF output from mutect2 with coverage and MQ annotation according to #64
* variant ids are now updated from simple VCF file
  ### Fixed
* Fixed a bug with sambamba depth coverage reporting wrong exon and panel coverage see #68
* The json output is now properly formatted using yapf
* Strelka rule doesn't filter out PASS variants anymore fixes issue #63

[1.12.0] - 2018-07-06
---------------------

Coverage report patch

Added
^^^^^


* Added a new script to retrieve coverage report for a list of gene(s) and transcripts(s)
* Added sambamba exon depth rule for coverage report
* Added a new entry in reference json for exon bed file, this file generated using: https://github.com/hassanfa/GFFtoolkit
  ### Changed
* sambamba_depth rule changed to sambama_panel_depth
* sambamba depth now has fix-mate-overlaps parameter enabled
* sambamba string filter changed to ``unmapped or mate\_is\_unmapped) and not duplicate and not failed\_quality\_control``.
* sambamba depth for both panel and exon work on picard flag (rmdup or mrkdup).
  ### Fixed
* Fixed sambamba panel depth rule for redundant coverage parameter

[1.11.0] - 2018-07-05
---------------------

create config patch for single and paired mode

Changed
^^^^^^^


* create_config is now accepting a paired|single mode instead of analysis json template (see help for changes). It is
  not backward compatible
  ### Added
* analysis_{paired single}.json for creating config. Analysis.json is now obsolete.
  ### Fixed
* A bug with writing output for analysis config, and creating the path if it doesn't exist.
* A bug with manta rule to correctly set output files in config.
* A bug that strelka was still included in sample analysis.

[1.10.0] - 2018-06-07
---------------------

Added
^^^^^


* Markduplicate flag to analysis config

[1.9.0] - 2018-06-04
--------------------

Added
^^^^^


* Single mode for vardict, manta, and mutect.
* merge type for tumor only
  ### Changed
* Single mode variant calling now has all variant calling rules
  ### Fixed
* run_analaysis now accepts workflows for testing pyrposes

[1.8.0] - 2018-06-01
--------------------

Changed
^^^^^^^


* picard create bed interval rule moved into collect hsmetric
* split bed is dependent on bam merge rule
* vardict env now has specific build rather than URL download (conda doesn't support URLs anymore)
  ### Fixed
* new logs and scripts dirs are not re-created if they are empty

[1.7.0] - 2018-05-31
--------------------

Added
^^^^^


* A source altered picard to generated more quality metrics output is added to installation and rules

[1.6.0] - 2018-05-30
--------------------

Added
^^^^^


* report subcommand for generating a pdf report from a json input file
* Added fastqc after removing adapter
  ### Changed
* Markduplicate now has both REMOVE and MARK (rmdup vs mrkdup)
* CollectHSMetrics now has more steps on PCT_TARGET_BASES

[1.5.0] - 2018-05-28
--------------------

Changed
^^^^^^^


* New log and script directories are now created for each re-run
  ### Fixed
* Picardtools' memory issue addressed for large samples

[1.4.0] - 2018-05-18
--------------------

Added
^^^^^


* single sample analysis mode
* alignment and insert size metrics are added to the workflow
  ### Changed
* collectqc and contest have their own rule for paired (tumor vs normal) and single (tumor only) sample.

[1.3.0] - 2018-05-13
--------------------

Added
^^^^^


* bed file for panel analysis is now mandatory to create analaysis config

[1.2.3] - 2018-05-13
--------------------

Changed
^^^^^^^


* vep execution path
* working directory for snakemake

[1.2.2] - 2018-05-04
--------------------

Added
^^^^^


* sbatch submitter and cluster config now has an mail field
  ### Changed
* ``create_config`` now only requires sample and output json. The rest are optional

[1.2.0] - 2018-05-02
--------------------

Added
^^^^^


* snakefile and cluster config in run analysis are now optional with a default value

[1.1.2] - 2018-04-27
--------------------

Fixed
^^^^^


* vardict installation was failing without conda-forge channel
* gatk installation was failing without correct jar file

[1.1.1] - 2018-04-27
--------------------

Fixed
^^^^^


* gatk-register tmp directory

[1.1.0] - 2018-04-26
--------------------

Added
^^^^^


* create config sub command added as a new feature to create input config file
* templates to generate a config file for analysis added
* code style template for YAPF input created. see: https://github.com/google/yapf
* vt conda env added

Changed
^^^^^^^


* install script changed to create an output config
* README updated with usage

Fixed
^^^^^


* fastq location for analysis config is now fixed
* lambda rules removed from cutadapt and fastq

[1.0.3-rc2] - 2018-04-18
------------------------

Added
^^^^^


* Added sbatch submitter to handle it outside snakemake
  ### Changed
* sample config file structure changed
* coding styles updated

[1.0.2-rc2] - 2018-04-17
------------------------

Added
^^^^^


* Added vt environment
  ### Fixed
* conda envs are now have D prefix instead of P (develop vs production)
* install_conda subcommand now accepts a proper conda prefix

[1.0.1-rc2] - 2018-04-16
------------------------

Fixed
^^^^^


* snakemake rules are now externally linked

[1.0.0-rc2] - 2018-04-16
------------------------

Added
^^^^^


* run_analysis subcommand
* Mutational Signature R script with CLI
* unittest to install_conda
* a method to semi-dynamically retrieve suitable conda env for each rule

Fixed
^^^^^


* install.sh updated with gatk and proper log output
* conda environments updated
* vardict now has its own environment and it should not raise anymore errors

[1.0.0-rc1] - 2018-04-05
------------------------

Added
^^^^^


* install.sh to install balsamic
* balsamic barebone cli
* subcommand to install required environments
* README.md updated with basic installation instructions

Fixed
^^^^^


* conda environment yaml files
