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
COPY varcall_py36.yaml ./varcall_py36.yaml

RUN conda env create --file varcall_py36.yaml -n varcall_py36 && \
    source activate varcall_py36 && \
    ln -s /usr/local/miniconda/envs/varcall_py36/lib/libreadline.so.7.0 /usr/local/miniconda/envs/varcall_py36/lib/libreadline.so.6 && \
    ln -s /usr/local/miniconda/envs/varcall_py36/lib/libreadline.so.7.0 /usr/local/miniconda/envs/varcall_py36/lib/libreadline.so.6.0 && \
    conda clean --index-cache --lock --tarballs -y
