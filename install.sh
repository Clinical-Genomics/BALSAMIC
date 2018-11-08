#!/bin/bash

B_RED='\033[0;31m';
B_GRN='\033[0;32m';
B_YLW='\033[1;33m';
B_NOCOL='\033[0m';
CONDAPREFIX=D
CONDADATE=$(date +%y%m%d)
CONDAPATH=/mnt/hds/proj/bioinfo/SERVER/miniconda/envs

#if [ $# -eq 0 ]; then
#  echo $"
#USAGE: $0 [-s CONDAPREFIX -d CONDADATE -p CONDAPATH]
#  1. Conda naming convention: [P,D]_[ENVNAME]_%DATE. P: Production, D: Development
#  2. Conda environment prefix: Path to conda env. e.g. /home/user/conda_env/
#  
#  -s CONDAPREFIX  Conda env name prefix. This will be P or D in the help above. 
#  -d CONDADATE    Conda env name suffix. This will be a suffix, by default it will be current date: yymmdd 
#  -p CONDAPATH    Conda env path prefix. See point 2 in help above.
#" >&2
#  exit 0
#fi

while getopts ":s:p:d:h" opt; do
  case $opt in
    s) sFlag=true;CONDAPREFIX=${OPTARG};;
    d) dFlag=true;CONDADATE=${OPTARG};;
    p) pFlag=true;CONDAPATH=${OPTARG};;
    h)
      echo $"
USAGE: $0 [-s CONDAPREFIX -d CONDADATE -p CONDAPATH]
  1. Conda naming convention: [P,D]_[ENVNAME]_%DATE. P: Production, D: Development
  2. Conda environment prefix: Path to conda env. e.g. /home/user/conda_env/
  
  -s CONDAPREFIX  Conda env name prefix. This will be P or D in the help above. 
  -d CONDADATE    Conda env name suffix. This will be a suffix, by default it will be current date: yymmdd 
  -p CONDAPATH    Conda env path prefix. See point 2 in help above.
" >&2
      exit 0
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires a argument." >&2
      exit 1
      ;;
  esac
done

# Check if conda exists
command -v conda > /dev/null 2>&1 || \
  { echo -e "${B_RED}conda command was not found. Please make sure conda is installed and it is in path. Aborting." >&2;\
  exit 1;
  }

# Conda env found
# Conda env naming convention: [P,D]_BALSAMIC_%DATE
# P: Production, D: Development 
env_name_suffix=${CONDADATE} 
env_name=${CONDAPREFIX}_BALSAMIC_${env_name_suffix}
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
balsamic install -s ${env_name_suffix} \
  --overwrite-env \
  --input-conda-yaml BALSAMIC/conda_yaml/D_Cancer-vardict.yaml \
  --input-conda-yaml BALSAMIC/conda_yaml/D_Cancer-Core.yaml \
  --input-conda-yaml BALSAMIC/conda_yaml/D_Cancer-py36.yaml \
  --input-conda-yaml BALSAMIC/conda_yaml/D_Cancer-py27.yaml \
  --input-conda-yaml BALSAMIC/conda_yaml/D_Cancer-vt.yaml \
  --env-dir-prefix ${CONDAPATH} \
  --packages-output-yaml ${BALSAMIC_ENVS}

gatk_env=`python -c 'from BALSAMIC.tools import get_conda_env; print(get_conda_env("BALSAMIC_env.yaml", "gatk"))'`

source activate ${gatk_env}

gatk-register BALSAMIC/assests/GenomeAnalysisTK.jar

echo -ne "${B_NOCOL}"
echo -e "\n${B_GRN}Copying custom Picard to relevant conda environment"
source activate ${env_name}
picard_PATH=BALSAMIC/assests/picard-2.18.5-4-g63d860e-SNAPSHOT-all.jar
picard_conda_env=`python -c 'from BALSAMIC.tools import get_conda_env; print(get_conda_env("BALSAMIC_env.yaml", "picard"))'`
picard_destination=`conda env export -n ${picard_conda_env} | grep prefix | cut -d" " -f 2`
cp $picard_PATH ${picard_destination}/share/
ln -s ${picard_destination}/share/picard-2.18.5-4-g63d860e-SNAPSHOT-all.jar  ${picard_destination}/share/picard-2.18.5.jar

echo -e "\n${B_YLW}Installing stargazer"
stargazer_path=BALSAMIC/assests/stargazer_5.2.2.tar.gz
R CMD INSTALL $stargazer_path

echo -e "\n${B_GRN}Install finished. Make sure you set reference.json and cluster.json."
echo -e "\n${B_GRN}To start working with BALSAMIC, run: source activate ${env_name}."

unset B_RED
unset B_GRN
unset B_YLW
unset B_NOCOL
