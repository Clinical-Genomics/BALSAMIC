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
COPY varcall_py27.yaml ./varcall_py27.yaml

RUN conda env create --file varcall_py27.yaml -n varcall_py27 \
    && conda clean --index-cache --lock --tarballs -y
