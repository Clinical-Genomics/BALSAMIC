#!/bin/bash
set -eo pipefail
shopt -s expand_aliases

_log="[$(date) $(whoami)] "
_red=${_log}'\033[0;31m';
_green=${_log}'\033[0;32m';
_yellow=${_log}'\033[1;33m';
_nocol='\033[0m';
_condaprefix=D

while getopts ":s:v:e:p:ch" opt; do
  case $opt in
    s) sFlag=true;_condaprefix=${OPTARG};;
    v) vFlag=true;_balsamic_ver=${OPTARG};;
    e) eFlag=true;_envsuffix=${OPTARG};;
    p) pFlag=true;_condapath=${OPTARG};;
    c) cFlag=true;;
    h)
      echo $"
USAGE: $0 [-s _condaprefix -v _balsamic_ver -p _condapath -c]
  1. Conda naming convention: [P,D,S]_[ENVNAME]_%DATE. P: Production, D: Development, S: Stage
  2. Conda environment prefix: Path to conda env. e.g. /home/user/conda_env/
  
  -s _condaprefix     Conda env name prefix. This will be P or D in the help above
  -v _balsamic_ver    Balsamic version tag to install (4.0.0+), or it could be the branch name
  -e _envsuffix       Balsamic conda env suffix. This will be added to the conda env name
  -p _condapath       Conda env path prefix. See point 2 in help above
  -c                  If set it will use Singularity container for conda instead 
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

if [[ -z $_balsamic_ver ]]
then
  echo -e "\n${_yellow}WARNING: No version or branch is set, master branch will be used.${_nocol}"
  _balsamic_ver="master"
fi

if [[ -z $_envsuffix ]]
then
  echo -e "\n${_yellow}WARNING: No conda env suffix provdided.${_nocol}"
fi

# Check if container flag is specified
if [[ $cFlag ]]
then
  current_conda_sif=${PWD}'/BALSAMIC/containers/BALSAMIC_miniconda3_4_6_14.sif'
  if [[ -f ${current_conda_sif} ]];
  then
    echo -e "\n${_green}Container for miniconda3 4.6.14 exists: ${current_conda_sif} ${_nocol}"
  else
    echo -e "\n${_green}Pulling a miniconda3 4.6.14 from shub://Clinical-Genomics/BALSAMIC:miniconda3_4_6_14.${_nocol}"
    singularity pull ${PWD}'/BALSAMIC/containers/BALSAMIC_miniconda3_4_6_14.sif' docker://hassanf/miniconda3_4.6.14 
  fi
  function conda() {
    singularity run --bind ${_condapath} ${PWD}'/BALSAMIC/containers/BALSAMIC_miniconda3_4_6_14.sif' conda "$@"
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

# Create conda environment
_env_name=${_condaprefix}_BALSAMIC${_envsuffix}
_balsamic_envs=${PWD}'/BALSAMIC_env.yaml'
_balsamic_ruledir=${PWD}'/BALSAMIC/'

echo -e "${_green}Creating conda env ${_env_name}${_nocol}"
conda env create --file BALSAMIC/conda/BALSAMIC-base.yaml --quiet --prefix ${_condapath}/${_env_name} --force

echo -e "${_green}Activating ${_env_name}${_nocol}"
source activate ${_env_name}

echo -e "${_green}Installing BALSAMIC from origin.${_nocol}"
pip install -U git+https://github.com/Clinical-Genomics/BALSAMIC@${_balsamic_ver}

# Pull Docker container
if [[ ${_balsamic_ver} =~ ^[0-9]+\.[0-9]+\.[0-9]+ ]]
then
  echo -e "${_green}Pulling version ${_balsamic_ver} of container.${_nocol}"
  container_version=release_v${_balsamic_ver} 
else
  echo -e "${_green}Pulling latest version of container.${_nocol}"
  container_version=latest
fi


_container_path=${PWD}"/BALSAMIC/containers/BALSAMIC_${container_version}.sif"
_docker_path=docker://hassanf/balsamic:${container_version}
echo -e "${_green}Grabbing container: ${_container_path} ${_nocol}"
singularity pull --force ${_container_path} ${_docker_path}

echo -e "\n${_green}Install finished. To start working with BALSAMIC, run: source activate ${_env_name}.${_nocol}"

unset _red
unset _green
unset _yellow
unset _nocol
