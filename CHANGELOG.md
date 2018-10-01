# Change Log

This change log will document the notable changes to this project in this file and it is following [Semantic
Versioning](https://semver.org/)

## [2.3.1] - 2018-09-25
### Added
- CollectHSmetric now run twice for before and after markduplicate

## [2.3.0] - 2018-09-25
### Changed
- Sample config file now includes a list of chromosomes in the panel bed file

### Fixed
- Non-matching chrom won't break the splitbed rule anymore
- collectqc rules now properly parse tab delimited metric files

## [2.2.0] - 2018-09-11
### Added
- Coverage plot to report
- target coverage file to report json
- post-cutadapt fastqc to collectqc
- A header to report pdf
- list of bioinfo tools used in the analysis added to report
### Changed
- VariantRep.R now accepts multiple inputs for each parameter (see help)
- AF values for MSKIMPACT config
### Fixed
- Output figure for coverageplot is now fully square :-)

## [2.1.0] - 2018-09-11
### Added
- normalized coverage plot script
- fastq file IO check for config creation
- added qos option to `balsamic run`
### Fixed
- Sambamba depth coverage parameters
- bug with picard markduplicate flag

## [2.0.2] - 2018-09-11
### Added
- Added qos option for setting qos to run jobs with a default value of low

## [2.0.1] - 2018-09-10
### Fixed
- Fixed package dependencies with vep and installation

## [2.0.0] - 2018-09-05
Variant reporter patch and cli update
### Added
- Added `balsamic config sample` and `balsamic config report` to generate run analysis and reporting config
- Added `VariantRep.R` script to information from merged variant table: variant summry, TMB, and much more
- Added a workflow for single sample mode alignment and QC only
- Added QC skimming script to qccollect to generate nicely formatted information from picard
### Changed
- Change to CLI for running and creating config
- Major overhaul to coverage report script. It's now simpler and more readable!
### Fixed
- Fixed sambamba depth to include mapping quality
- Markduplicate now is now by default on marking mode, and will NOT remove duplicates
- Minor formatting and script beautification happened

## [1.13.1] - 2018-08-17
### Fixed
- fixed a typo in MSKMVL config
- fixed a bug in strelka\_simple for correct column orders

## [1.13.0] - 2018-08-10
### Added
- rule for all three variant callers for paired analysis now generate a simple VCF file
- rule for all three variant callers for paired analysis to convert VCF into table format
- MVL config file and MVL annotation to VCF calls for SNV/INDEL callers
- CALLER annotation added to SNV/INDEL callers
- exome specific option for strelka paired
- create\_config subcommand is now more granular, it accepts all enteries from sample.json as commandline arguments
- Added tabQuery to the assests as a tool to query the tabulated output of summarized VCF
- Added MQ annotation field to Mutect2 output see #67
### Changed
- Leaner VCF output from mutect2 with coverage and MQ annotation according to #64
- variant ids are now updated from simple VCF file
### Fixed
- Fixed a bug with sambamba depth coverage reporting wrong exon and panel coverage see #68
- The json output is now properly formatted using yapf
- Strelka rule doesn't filter out PASS variants anymore fixes issue #63

## [1.12.0] - 2018-07-06
Coverage report patch
### Added
- Added a new script to retrieve coverage report for a list of gene(s) and transcripts(s)
- Added sambamba exon depth rule for coverage report
- Added a new entry in reference json for exon bed file, this file generated using: https://github.com/hassanfa/GFFtoolkit
### Changed
- sambamba\_depth rule changed to sambama\_panel\_depth
- sambamba depth now has fix-mate-overlaps parameter enabled
- sambamba string filter changed to `unmapped or mate\_is\_unmapped) and not duplicate and not failed\_quality\_control`.
- sambamba depth for both panel and exon work on picard flag (rmdup or mrkdup).
### Fixed
- Fixed sambamba panel depth rule for redundant coverage parameter

## [1.11.0] - 2018-07-05
create config patch for single and paired mode
### Changed
- create\_config is now accepting a paired|single mode instead of analysis json template (see help for changes). It is
  not backward compatible
### Added
- analysis\_{paired single}.json for creating config. Analysis.json is now obsolete.
### Fixed
- A bug with writing output for analysis config, and creating the path if it doesn't exist.
- A bug with manta rule to correctly set output files in config.
- A bug that strelka was still included in sample analysis.

## [1.10.0] - 2018-06-07
### Added
- Markduplicate flag to analysis config

## [1.9.0] - 2018-06-04
### Added
- Single mode for vardict, manta, and mutect.
- merge type for tumor only
### Changed
- Single mode variant calling now has all variant calling rules
### Fixed
- run\_analaysis now accepts workflows for testing pyrposes

## [1.8.0] - 2018-06-01
### Changed
- picard create bed interval rule moved into collect hsmetric
- split bed is dependent on bam merge rule
- vardict env now has specific build rather than URL download (conda doesn't support URLs anymore)
### Fixed
- new logs and scripts dirs are not re-created if they are empty

## [1.7.0] - 2018-05-31
### Added
- A source altered picard to generated more quality metrics output is added to installation and rules

## [1.6.0] - 2018-05-30
### Added
- report subcommand for generating a pdf report from a json input file
- Added fastqc after removing adapter
### Changed
- Markduplicate now has both REMOVE and MARK (rmdup vs mrkdup)
- CollectHSMetrics now has more steps on PCT\_TARGET\_BASES

## [1.5.0] - 2018-05-28
### Changed
- New log and script directories are now created for each re-run
### Fixed
- Picardtools' memory issue addressed for large samples

## [1.4.0] - 2018-05-18
### Added
- single sample analysis mode
- alignment and insert size metrics are added to the workflow
### Changed
- collectqc and contest have their own rule for paired (tumor vs normal) and single (tumor only) sample.

## [1.3.0] - 2018-05-13
### Added
- bed file for panel analysis is now mandatory to create analaysis config

## [1.2.3] - 2018-05-13
### Changed
- vep execution path
- working directory for snakemake

## [1.2.2] - 2018-05-04
### Added
- sbatch submitter and cluster config now has an mail field
### Changed
- `create_config` now only requires sample and output json. The rest are optional

## [1.2.0] - 2018-05-02
### Added
- snakefile and cluster config in run analysis are now optional with a default value

## [1.1.2] - 2018-04-27
### Fixed
- vardict installation was failing without conda-forge channel
- gatk installation was failing without correct jar file

## [1.1.1] - 2018-04-27
### Fixed
- gatk-register tmp directory

## [1.1.0] - 2018-04-26
### Added
- create config sub command added as a new feature to create input config file
- templates to generate a config file for analysis added
- code style template for YAPF input created. see: https://github.com/google/yapf
- vt conda env added

### Changed
- install script changed to create an output config
- README updated with usage

### Fixed
- fastq location for analysis config is now fixed
- lambda rules removed from cutadapt and fastq

## [1.0.3-rc2] - 2018-04-18
### Added
- Added sbatch submitter to handle it outside snakemake
### Changed
- sample config file structure changed
- coding styles updated

## [1.0.2-rc2] - 2018-04-17
### Added
- Added vt environment
### Fixed
- conda envs are now have D prefix instead of P (develop vs production)
- install_conda subcommand now accepts a proper conda prefix

## [1.0.1-rc2] - 2018-04-16
### Fixed
- snakemake rules are now externally linked

## [1.0.0-rc2] - 2018-04-16
### Added
- run_analysis subcommand
- Mutational Signature R script with CLI
- unittest to install_conda
- a method to semi-dynamically retrieve suitable conda env for each rule

### Fixed
- install.sh updated with gatk and proper log output
- conda environments updated
- vardict now has its own environment and it should not raise anymore errors

## [1.0.0-rc1] - 2018-04-05
### Added
- install.sh to install balsamic
- balsamic barebone cli
- subcommand to install required environments
- README.md updated with basic installation instructions

### Fixed
- conda environment yaml files
