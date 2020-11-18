FROM continuumio/miniconda3:4.8.2-alpine 

LABEL base_image="continuumio/miniconda3:4.8.2-alpine"
LABEL about.home="https://github.com/Clinical-Genomics/vogue"
LABEL about.documentation="https://vogue.readthedocs.io/"
LABEL about.license="MIT License (MIT)"
LABEL about.maintainer="Hassan Foroughi hassan dot foroughi at scilifelab dot se" 
LABEL about.description="Bioinformatic analysis pipeline for somatic mutations in cancer"
LABEL about.version="6.0.1"

ENV PATH="/opt/conda/bin/:${PATH}"

# Copy all project files
COPY align.yaml ./align.yaml

RUN conda env create -n align_qc --file align.yaml && \
    source activate align_qc && \
    picard_jar=picard-2.23.2-201-g922891d-SNAPSHOT-all.jar && \
    picard_PATH=BALSAMIC/assets/${picard_jar} && \
    picard_destination=/usr/local/miniconda/envs/align_qc/share/ && \
    cp $picard_PATH ${picard_destination} && \
    ln -s ${picard_destination}/${picard_jar} ${picard_destination}/picard.jar; \
    ln -s /usr/local/miniconda/envs/align_qc/lib/libreadline.so.7.0 /usr/local/miniconda/envs/align_qc/lib/libreadline.so.6 && \
    ln -s /usr/local/miniconda/envs/align_qc/lib/libreadline.so.7.0 /usr/local/miniconda/envs/align_qc/lib/libreadline.so.6.0
    && conda clean --index-cache --lock --tarballs -y
