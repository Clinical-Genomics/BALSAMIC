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
- Activate BALSAMIC conda environment and run setup.py to install BALSAMIC itself


