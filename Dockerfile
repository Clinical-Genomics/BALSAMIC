FROM continuumio/miniconda3:4.9.2-alpine

LABEL base_image="continuumio/miniconda3:4.9.2-alpine"
LABEL about.home="https://github.com/Clinical-Genomics/BALSAMIC"
LABEL about.documentation="https://balsamic.readthedocs.io/"
LABEL about.license="MIT License (MIT)"
LABEL about.maintainer="Hassan Foroughi hassan dot foroughi at scilifelab dot se" 
LABEL about.description="Bioinformatic analysis pipeline for somatic mutations in cancer"
LABEL about.version="7.2.6"

ENV PATH="/opt/conda/bin/:${PATH}"
RUN apk add --no-cache bash gcc git libc-dev zlib-dev linux-headers

ARG CONTAINER_NAME
ARG WORK_DIR=project

# Copy all project files
COPY . /${WORK_DIR}

RUN cd /${WORK_DIR}/BALSAMIC/containers/${CONTAINER_NAME}/ && /bin/sh ${CONTAINER_NAME}.sh ${CONTAINER_NAME}
RUN if [ "$CONTAINER_NAME" = "balsamic" ]; then cd /project && pip install --no-cache-dir . ; fi

# Clean work environment
RUN rm -rf /project && conda clean --all --yes
