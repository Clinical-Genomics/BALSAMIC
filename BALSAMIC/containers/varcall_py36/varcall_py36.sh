conda env create -n ${1} --file ${1}.yaml
source activate ${1}
ENV_PATH=/opt/conda/envs
ln -s ${ENV_PATH}/${1}/lib/libreadline.so.7.0 ${ENV_PATH}/${1}/lib/libreadline.so.6
ln -s ${ENV_PATH}/${1}/lib/libreadline.so.7.0 ${ENV_PATH}/${1}/lib/libreadline.so.6.0
