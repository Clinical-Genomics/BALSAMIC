#!/bin/bash
set -eo pipefail
shopt -s expand_aliases

_log="[$(date) $(whoami)] "
_red=${_log}'\033[0;31m';
_green=${_log}'\033[0;32m';
_yellow=${_log}'\033[1;33m';
_nocol='\033[0m';

function usage() {
  echo $"
USAGE: [ -a <T|TN> -c _condaenv -m <run|config|all> -t <panel|WGS> -r ]

  -w Analysis workflow: balsamic [default option], balsamic-umi or balsamic-qc
  -a [required] T: Tumor only, TN: tumor normal
  -c Conda environment where BALSAMIC is installed. If not specified, it will use current environment.
  -m [required] config: only create config file, run: create config file and start analysis
  -t [required] panel: target sequencing workflow (also includes WES), WGS: whole genome sequencing workflow
  -d Analysis dir path, if it doesn't exist it will be created
  -r Flag. Set to submit jobs instead of running in dry mode
  -h Show this help and exit
"
}

while getopts ":w:a:m:c:t:d:r" opt; do
  case ${opt} in
    w)
      _analysis_workflow=${OPTARG}
      echo "analysis workflow set to" "${OPTARG}"
      ;;
    a)
      _analysis=${OPTARG}
      echo "analysis set to" "${OPTARG}"
      [[ $_analysis == 'T' || $_analysis == 'TN' ]] || ( usage >&2; exit 1)
      ;;
    c)
      _condaenv=${OPTARG}
      echo "conda environment set to" "${OPTARG}"
      ;;
    m)
      _startmode=${OPTARG}
      echo "start mode set to" "${OPTARG}"
      [[ $_startmode == 'config' || $_startmode == 'run' || $_startmode == 'all' ]] || ( usage >&2; exit 1)
      ;;
    t)
      _ngstype=${OPTARG}
      echo "workflow set to " "${OPTARG}"
      [[ $_ngstype == 'panel' || $_ngstype == 'WGS' ]] || ( usage >&2; exit 1)
      ;;
    d)
      _analysis_dir=${OPTARG}
      echo "analysis dir set to " "${OPTARG}"
      ;;
    r)
      rFlag=true;
      ;;
    *) echo "Invalid option: -${OPTARG}" >&2; usage >&2; exit 1;;
  esac
done

if [[ ${_ngstype} == "panel" ]]; then
  _panel_option='-p tests/test_data/references/panel/panel.bed'
else
  _panel_option=''
fi

if [[ ! -z ${_condaenv} ]]; then
  source activate ${_condaenv}
fi

if [[ -z ${_analysis_dir} ]]; then
  _analysis_dir='run_tests/'
  echo "analysis dir set to " "${_analysis_dir}"
fi

# Make sure _analysis_dir exists
mkdir -p ${_analysis_dir}

_genome_ver=hg19
_workflow=balsamic
_cluster_config=BALSAMIC/config/cluster.json
_balsamic_cache=/home/proj/stage/cancer/balsamic_cache
_tumor_fastq=tests/test_data/fastq/S1_R_1.fastq.gz
_normal_fastq=tests/test_data/fastq/S2_R_1.fastq.gz
_analysis_config=${_analysis_dir}'/'${_analysis}_${_ngstype}'/'${_analysis}_${_ngstype}'.json'

if [[ ! -z ${_analysis_workflow} ]]; then
  _workflow=${_analysis_workflow}
fi

if [[ ! -z ${rFlag} ]]; then
  _run_analysis="-r"
fi

if [[ ${_analysis} == "TN" ]]; then
  _normal_option="-n ${_normal_fastq}"
else
  _normal_option=" "
fi

function balsamic_config() {
set -x
  balsamic --loglevel INFO config case \
    -w ${_workflow} \
    -t ${_tumor_fastq} \
    ${_normal_option} \
    --case-id ${_analysis}_${_ngstype} \
    --analysis-dir ${_analysis_dir} \
    ${_panel_option} \
    --balsamic-cache ${_balsamic_cache}
}

balsamic_run() {
  balsamic --loglevel INFO run analysis \
    -s ${_analysis_config} \
    -c ${_cluster_config} \
    --benchmark \
    --account development ${_run_analysis}
}

if [[ $_startmode == 'config' ]]; then
  balsamic_config
elif [[ $_startmode == 'run' ]]; then
  balsamic_run
else
  balsamic_config
  balsamic_run
fi
