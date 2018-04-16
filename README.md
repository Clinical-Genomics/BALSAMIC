# BALSAMIC
Bioinformatic Analysis pipeLine for SomAtic MutatIons in Cancer

## Requirments

### Software

BALSAMIC requires conda 4.3+ for automated installation. For detailed software and python requirments please see
```requirments.txt``` and ```BALSAMIC/install/conda_yaml/*yaml``` 

## Installation

BALSAMIC requires multiple conda environments. ```install.sh``` will automatically install BALSAMIC in the following
order:

- Create a python 3.6 conda environment following the conda environment naming convention from
  [https://github.com/Clinical-Genomics/development/blob/master/conda/conda_conventions.md](https://github.com/Clinical-Genomics/development/blob/master/conda/conda_conventions.md)
based on the config file: ```BALSAMIC/conda_yaml/BALSAMIC.yaml``` 
- Create conda environments required for BALSAMIC to run properly
- Install BALSAMIC
- Install gatk


## Usage

```rule_directory``` within ```config_sample.json``` has to be set correctly for VariantCalling snakefile to work
properly. If it is set correctly, the analysis can be initiated using the following command:

```bash
balsamic run_analysis -S VariantCalling -c  config_sample.json -s config_cluster.json
```

Where VariantCalling, config_sample.json, and config_cluster are snakefile, sample config file, and cluster config file
respectively.

The template for structure of config_sample.json and config_cluster.json can be found inside
```BALSAMIC/config_files```

