# BALSAMIC
Bioinformatic Analysis pipeLine for SomAtic MutatIons in Cancer

## Requirments

### Software

BALSAMIC requires conda 4.3+ for automated installation. For detailed software and python requirments please see
```requirments.txt``` and ```BALSAMIC/install/conda_yaml/*ysml``` 

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

For now snakemake rules are not linked externally, so variant calling has to run within ```workflows``` directory:

```bash
balsamic run_analysis -S VariantCalling -c  ../config_files/config_sample.json -s ../config_files/config_cluster.json
```

The template for structure of config_sample.json and config_cluster.json can be found inside
```BALSAMIC/config_files``can be found inside ```BALSAMIC/config_files```

 


