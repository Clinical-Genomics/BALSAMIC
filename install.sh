#!/bin/bash

B_RED='\033[0;31m';
B_GRN='\033[0;32m';
B_YLW='\033[1;33m';
B_NOCOL='\033[0m';

# Check if conda exists
command -v conda > /dev/null 2>&1 || \
  { echo -e "${B_RED}conda command was not found. Please make sure conda is installed and it is in path. Aborting." >&2;\
  exit 1;
  }

# Conda env found
# Conda env naming convention: [P,D]_BALSAMIC_%DATE
# P: Production, D: Development
env_name=P_BALSAMIC_$(date +%y%m%d)

echo -ne "${B_NOCOL}"
echo -e "\n${B_GREEN}Creating conda env ${env_name}"
echo -e "\n${B_YLW}\tconda env create -f BALSAMIC/conda_yaml/BALSAMIC.yaml -c bioconda -c conda-forge --quiet --name ${env_name} python=3.6"

echo -ne "${B_NOCOL}"
echo -e "\n${B_GREEN}Update .bashrc"
echo -e "\n${B_YLW}\techo 'export PYTHONPATH=\$PYTHONPATH:/home/hassan.foroughi/repo/BALSAMIC/' >> ~/.bashrc"

echo -ne "${B_NOCOL}"
echo -e "\n${B_GREEN}Activating ${env_name}"
echo -e "\n${B_YLW}\tsource activate ${env_name}"

echo -ne "${B_NOCOL}"
echo -e "\n${B_GREEN}Installing snakemake"
echo -e "\n${B_YLW}\tconda install -c bioconda snakemake=4.8.0"

echo -ne "${B_NOCOL}"
echo -e "\n${B_GREEN}Installing pip"
echo -e "\n${B_YLW}\tconda install -c bioconda pip=9.0.3"

echo -ne "${B_NOCOL}"
echo -e "\n${B_GREEN}Installing requirments.txt"
echo -e "\n${B_YLW}\tpip install requirments.txt"

echo -ne "${B_NOCOL}"
echo -e "\n${B_GREEN}Installing BALSAMIC"
echo -e "\n${B_YLW}\tpip install ../"

echo -ne "${B_NOCOL}"
echo -e "\n${B_GREEN}Registering GATK 3.8.0"
echo -e "\n${B_YLW}\tsource activate GATKENV"
echo -e "\n${B_YLW}\tgatk-register gatk_3.8.tar.bz2 "

unset B_RED
unset B_GRN
unset B_YLW
unset B_NOCOL
