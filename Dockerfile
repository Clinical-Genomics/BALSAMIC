FROM continuumio/miniconda3:4.9.2-alpine

LABEL base_image="continuumio/miniconda3:4.9.2-alpine"
LABEL about.home="https://github.com/Clinical-Genomics/BALSAMIC"
LABEL about.documentation="https://balsamic.readthedocs.io/"
LABEL about.license="MIT License (MIT)"
LABEL about.maintainer="Hassan Foroughi hassan dot foroughi at scilifelab dot se" 
LABEL about.description="Bioinformatic analysis pipeline for somatic mutations in cancer"
LABEL about.version="7.2.5"

ENV PATH="/opt/conda/bin/:${PATH}"
RUN apk add --no-cache bash gcc git 

ARG CONTAINER_NAME

# Copy all project files
COPY BALSAMIC/containers/${CONTAINER_NAME}/${CONTAINER_NAME}.yaml ./${CONTAINER_NAME}.yaml
COPY BALSAMIC/containers/${CONTAINER_NAME}/${CONTAINER_NAME}.sh ./${CONTAINER_NAME}.sh

RUN /bin/sh ${CONTAINER_NAME}.sh ${CONTAINER_NAME} && conda clean --all --yes 

