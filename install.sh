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
env_name_suffix=_180410 #_$(date +%y%m%d)
env_name=D_BALSAMIC${env_name_suffix}
BALSAMIC_ENVS=${PWD}'/BALSAMIC_env.yaml'
BALSAMIC_RULEDIR=${PWD}'/BALSAMIC/'

echo -ne "${B_NOCOL}"
echo -e "\n${B_GRN}Writing BALSAMIC/config/install.json ${env_name}"
cat > BALSAMIC/config/install.json << EOF
{
    "conda_env_yaml": "${BALSAMIC_ENVS}",
    "rule_directory": "${BALSAMIC_RULEDIR}"
}
EOF
echo -e "${B_YLW}"
cat BALSAMIC/config/install.json

echo -ne "${B_NOCOL}"
echo -e "\n${B_GRN}Creating conda env ${env_name}"
echo -e "\n${B_YLW}\tconda env create -f BALSAMIC/conda_yaml/BALSAMIC.yaml --quiet --name ${env_name} --force"
conda env create -f BALSAMIC/conda_yaml/BALSAMIC.yaml --quiet --prefix /mnt/hds/proj/bioinfo/SERVER/miniconda/envs/${env_name} --force

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
echo -e "\n${B_YLW}\tbalsamic install_env --packages-output-yaml ${BALSAMIC_ENVS} -s ${env_name_suffix} -i BALSAMIC/conda_yaml/D_Cancer-vardict.yaml -i BALSAMIC/conda_yaml/D_Cancer-Core.yaml -i BALSAMIC/conda_yaml/D_Cancer-py36.yaml -i BALSAMIC/conda_yaml/D_Cancer-py27.yaml -i BALSAMIC/conda_yaml/D_Cancer-vt.yaml -o"
echo -ne "${B_NOCOL}"
balsamic install_env -s ${env_name_suffix} \
  --overwrite-env \
  --input-conda-yaml BALSAMIC/conda_yaml/D_Cancer-vardict.yaml \
  --input-conda-yaml BALSAMIC/conda_yaml/D_Cancer-Core.yaml \
  --input-conda-yaml BALSAMIC/conda_yaml/D_Cancer-py36.yaml \
  --input-conda-yaml BALSAMIC/conda_yaml/D_Cancer-py27.yaml \
  --input-conda-yaml BALSAMIC/conda_yaml/D_Cancer-vt.yaml \
  --env-dir-prefix /mnt/hds/proj/bioinfo/SERVER/miniconda/envs \
  --packages-output-yaml ${BALSAMIC_ENVS}

gatk_env=`python -c 'from BALSAMIC.tools import get_conda_env; print(get_conda_env("BALSAMIC_env.yaml", "gatk"))'`

source activate ${gatk_env}

echo -e "\n${B_YLW}Updating gatk tmp directory to ${CONDA_PREFIX}/gatk_tmp"
sed -i 's/\/tmp\/gatk/\$\{CONDA_PREFIX\}\/gatk_tmp/g' ${CONDA_PREFIX}/opt/gatk-3.8/gatk-register.sh

gatk-register BALSAMIC/install/gatk_3.8.tar.bz2

unset B_RED
unset B_GRN
unset B_YLW
unset B_NOCOL
