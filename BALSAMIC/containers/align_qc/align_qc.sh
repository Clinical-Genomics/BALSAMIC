conda env create -n ${1} --file ${1}.yaml

#&& \
#    source activate align_qc && \
#    picard_jar=picard-2.23.2-201-g922891d-SNAPSHOT-all.jar && \
#    picard_PATH=BALSAMIC/assets/${picard_jar} && \
#    picard_destination=/usr/local/miniconda/envs/align_qc/share/ && \
#    cp $picard_PATH ${picard_destination} && \
#    ln -s ${picard_destination}/${picard_jar} ${picard_destination}/picard.jar; \
#    ln -s /usr/local/miniconda/envs/align_qc/lib/libreadline.so.7.0 /usr/local/miniconda/envs/align_qc/lib/libreadline.so.6 && \
#    ln -s /usr/local/miniconda/envs/align_qc/lib/libreadline.so.7.0 /usr/local/miniconda/envs/align_qc/lib/libreadline.so.6.0
#    && conda clean --index-cache --lock --tarballs -y
