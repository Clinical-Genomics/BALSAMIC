FROM continuumio/miniconda3:4.8.2-alpine 

LABEL base_image="continuumio/miniconda3:4.8.2-alpine"
LABEL about.home="https://github.com/Clinical-Genomics/BALSAMIC"
LABEL about.documentation="https://balsamic.readthedocs.io/"
LABEL about.license="MIT License (MIT)"
LABEL about.maintainer="Hassan Foroughi hassan dot foroughi at scilifelab dot se" 
LABEL about.description="Bioinformatic analysis pipeline for somatic mutations in cancer"
LABEL about.version="6.0.2"

ENV PATH="/opt/conda/bin/:${PATH}"

ARG CONTAINER_NAME

# Copy all project files
COPY BALSAMIC/containers/${CONTAINER_NAME}/${CONTAINER_NAME}.yaml ./${CONTAINER_NAME}.yaml
COPY BALSAMIC/containers/${CONTAINER_NAME}/${CONTAINER_NAME}.sh ./${CONTAINER_NAME}.sh

RUN /bin/sh ${CONTAINER_NAME}.sh ${CONTAINER_NAME} && conda clean --all --yes 
