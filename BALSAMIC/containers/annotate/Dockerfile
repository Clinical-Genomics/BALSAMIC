FROM continuumio/miniconda3:4.10.3-alpine

LABEL base.image="continuumio/miniconda3:4.10.3-alpine"
LABEL about.home="https://github.com/Clinical-Genomics/BALSAMIC"
LABEL about.documentation="https://balsamic.readthedocs.io/"
LABEL about.license="MIT License (MIT)"
LABEL about.maintainer="Hassan Foroughi hassan dot foroughi at scilifelab dot se" 
LABEL about.description="Bioinformatic analysis pipeline for somatic mutations in cancer"

ENV PATH="/opt/conda/bin/:${PATH}"

RUN apk add --no-cache bash gcc git zlib-dev musl-dev

ARG WORK_DIR=project
ARG CONTAINER_NAME

# Copy all project files
COPY . /${WORK_DIR}

RUN cd /${WORK_DIR}/BALSAMIC/containers/${CONTAINER_NAME}/ && /bin/sh ${CONTAINER_NAME}.sh ${CONTAINER_NAME}

# Clean work environment
RUN rm -rf /${WORK_DIR} && conda clean --all --yes
