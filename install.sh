#!/bin/bash
set -eo pipefail
shopt -s expand_aliases

_log="[$(date) $(whoami)] "
_red=${_log}'\033[0;31m';
_green=${_log}'\033[0;32m';
_yellow=${_log}'\033[1;33m';
_nocol='\033[0m';
_condaprefix=D
_condadate=$(date +%y%m%d)

#if [ $# -eq 0 ]; then
#  echo $"
#USAGE: $0 [-s _condaprefix -d _condadate -p _condapath -c]
#  1. Conda naming convention: [P,D]_[ENVNAME]_%DATE. P: Production, D: Development
#  2. Conda environment prefix: Path to conda env. e.g. /home/user/conda_env/
#  
#  -s _condaprefix  Conda env name prefix. This will be P or D in the help above. 
#  -d _condadate    Conda env name suffix. This will be a suffix, by default it will be current date: yymmdd 
#  -p _condapath    Conda env path prefix. See point 2 in help above.
#  -c If set it will use Singularity container for conda instead 
#" >&2
#  exit 0
#fi

while getopts ":s:p:d:ch" opt; do
  case $opt in
    s) sFlag=true;_condaprefix=${OPTARG};;
    d) dFlag=true;_condadate=${OPTARG};;
    p) pFlag=true;_condapath=${OPTARG};;
    c) cFlag=true;;
    h)
      echo $"
USAGE: $0 [-s _condaprefix -d _condadate -p _condapath -c]
  1. Conda naming convention: [P,D]_[ENVNAME]_%DATE. P: Production, D: Development
  2. Conda environment prefix: Path to conda env. e.g. /home/user/conda_env/
  
  -s _condaprefix  Conda env name prefix. This will be P or D in the help above. 
  -d _condadate    Conda env name suffix. This will be a suffix, by default it will be current date: yymmdd 
  -p _condapath    Conda env path prefix. See point 2 in help above.
  -c If set it will use Singularity container for conda instead 
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

if [[ -z $_condapath  ]]
then
  echo -e "\n${_red}No conda env path provided. Exiting!${_nocol}"
  exit 1
fi

# Check if container flag is specified
if [[ $cFlag ]]
then
  echo -e "\n${_green}Pulling a miniconda3 4.6.14 from shub://Clinical-Genomics/BALSAMIC:miniconda3_4_6_14.${_nocol}"
  singularity pull shub://Clinical-Genomics/BALSAMIC:miniconda3_4_6_14
  function conda() {
    singularity run --bind ${_condapath} BALSAMIC_miniconda3_4_6_14.sif conda "$@"
  }
fi

# Check if conda exists
if [[ -z $cFlag ]]
then
  command -v conda > /dev/null 2>&1 || \
    { >&2 echo -e "${_red}conda command was not found. Please make sure conda is installed and it is in path. Aborting.";\
      >&2 echo -e "${_red}If you want to installed without conda command available, consider using -c flag.";\
    exit 1;
    }
fi


# Conda env found
# Conda env naming convention: [P,D]_BALSAMIC_%DATE
# P: Production, D: Development 
_env_name_suffix=_${_condadate} 
_env_name=${_condaprefix}_BALSAMIC-base${_env_name_suffix}
_balsamic_envs=${PWD}'/BALSAMIC_env.yaml'
_balsamic_ruledir=${PWD}'/BALSAMIC/'

echo -e "${_green}Writing BALSAMIC/config/install.json ${_env_name}${_nocol}"
cat > BALSAMIC/config/install.json << EOF
{
    "conda_env_yaml": "${_balsamic_envs}",
    "rule_directory": "${_balsamic_ruledir}"
}
EOF
cat BALSAMIC/config/install.json

echo -e "${_green}Creating conda env ${_env_name}${_nocol}"
conda env create -f BALSAMIC/conda_yaml/BALSAMIC-base.yaml --quiet --prefix ${_condapath}/${_env_name} --force

echo -e "${_green}Activating ${_env_name}${_nocol}"
echo $_env_name
source activate ${_env_name}

echo -e "${_green}Installing BALSAMIC${_nocol}"
echo -e "${_yellow}\tpip install --editable .${_nocol}"
pip install -r requirements.txt --editable .

echo -e "${_green}Installting environments for the workflow.${_nocol}"

balsamic install -s ${_env_name_suffix} \
  --overwrite-env \
  --env-type ${_condaprefix} \
  --input-conda-yaml BALSAMIC/conda_yaml/BALSAMIC-py27.yaml \
  --input-conda-yaml BALSAMIC/conda_yaml/BALSAMIC-py36.yaml \
  --env-dir-prefix ${_condapath} \
  --packages-output-yaml ${_balsamic_envs}

gatk_env=`python -c 'from BALSAMIC.tools import get_conda_env; print(get_conda_env("BALSAMIC_env.yaml", "gatk"))'`

source activate ${gatk_env}
gatk3-register BALSAMIC/assets/GenomeAnalysisTK.jar

echo -e "${_green}Copying custom Picard to relevant conda environment.${_nocol}"
source activate ${_env_name}
picard_PATH=BALSAMIC/assets/picard-2.18.11-3-gc6e797f-SNAPSHOT-all.jar
picard_conda_env=`python -c 'from BALSAMIC.tools import get_conda_env; print(get_conda_env("BALSAMIC_env.yaml", "picard"))'`
picard_destination=${picard_conda_env}/share/
cp $picard_PATH ${picard_destination}
# link picard from assets to conda's share path
ln -s ${picard_destination}/picard-2.18.11-3-gc6e797f-SNAPSHOT-all.jar  ${picard_destination}/picard-2.18.11.jar

echo -e "\n${_green}Install finished. Make sure you set reference.json and cluster.json.${_nocol}"
echo -e "\n${_green}To start working with BALSAMIC, run: source activate ${_env_name}.${_nocol}"

unset _red
unset _green
unset _yellow
unset _nocol
