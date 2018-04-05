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
echo -e "\n${B_GRN}Creating conda env ${env_name}"
echo -e "\n${B_YLW}\tconda env create -f BALSAMIC/conda_yaml/BALSAMIC.yaml --quiet --name ${env_name} --force"
conda env create -f BALSAMIC/conda_yaml/BALSAMIC.yaml --quiet --name ${env_name} --force

echo -ne "${B_NOCOL}"
echo -e "\n${B_GRN}Activating ${env_name}"
echo -e "\n${B_YLW}\tsource activate ${env_name}"
source activate ${env_name}

echo -ne "${B_NOCOL}"
echo -e "\n${B_GRN}Installing requirments.txt"
echo -e "\n${B_YLW}\tpip install -r requirments.txt"
pip install -r requirments.txt

echo -ne "${B_NOCOL}"
echo -e "\n${B_GRN}Installing BALSAMIC"
echo -e "\n${B_YLW}\tpip install --editable ."
pip install --editable .

echo -ne "${B_NOCOL}"
echo -e "\n${B_GRN}Installting environments for the workflow"
echo -e "\n${B_YLW}\tbalsamic install_conda -i BALSAMIC/conda_yaml/P_Cancer-Core.yaml -s _$(date +%y%m%d)  -i BALSAMIC/conda_yaml/P_Cancer-py36.yaml -i BALSAMIC/conda_yaml/P_Cancer-py27.yaml -o"
echo -ne "${B_NOCOL}"
balsamic install_conda -s _$(date +%y%m%d) \
  --input-conda-yaml BALSAMIC/conda_yaml/P_Cancer-Core.yaml \
  --input-conda-yaml BALSAMIC/conda_yaml/P_Cancer-py36.yaml \
  --input-conda-yaml BALSAMIC/conda_yaml/P_Cancer-py27.yaml \
  --overwrite-env

unset B_RED
unset B_GRN
unset B_YLW
unset B_NOCOL
