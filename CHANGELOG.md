# Change Log

This change log will document the notable changes to this project in this file and it is following [Semantic
Versioning](https://semver.org/)

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
- create_config now only requires sample and output json. The rest are optional

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
