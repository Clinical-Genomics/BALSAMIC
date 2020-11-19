conda env create -n ${1} --file ${1}.yaml && \
    source activate ${1} && \
    ln -s /usr/local/miniconda/envs/${1}/lib/libreadline.so.7.0 /usr/local/miniconda/envs/${1}/lib/libreadline.so.6 && \
    ln -s /usr/local/miniconda/envs/${1}/lib/libreadline.so.7.0 /usr/local/miniconda/envs/${1}/lib/libreadline.so.6.0
