# Change Log

This change log will document the notable changes to this project in this file and it is following [Semantic
Versioning](https://semver.org/)

## [1.0.0-rc2] - 2018-04-05
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
