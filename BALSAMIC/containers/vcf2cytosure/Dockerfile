FROM continuumio/miniconda3:4.10.3-alpine

LABEL base.image="continuumio/miniconda3:4.10.3-alpine"
LABEL about.home="https://github.com/Clinical-Genomics/BALSAMIC"
LABEL about.documentation="https://balsamic.readthedocs.io/"
LABEL about.license="MIT License (MIT)"
LABEL about.maintainer="Ashwini Jeggari ashwini dot jeggari at scilifelab dot se"
LABEL about.description="Bioinformatic analysis pipeline for somatic mutations in cancer"

ARG CONTAINER_NAME
ARG WORK_DIR=project
ENV PATH="/opt/${CONTAINER_NAME}/bin:${PATH}"
ENV PATH="/opt/conda/bin/:${PATH}"
ENV PYTHONPATH="/opt/${CONTAINER_NAME}"

COPY . /${WORK_DIR}

RUN apk add --no-cache bash gcc git python3 gzip

RUN cd /opt \
    && git clone https://github.com/NBISweden/vcf2cytosure.git \
    && cd /opt/${CONTAINER_NAME}/  \
    && git checkout v0.8 \
    && pip install --no-cache-dir .

RUN cd /${WORK_DIR}/BALSAMIC/containers/${CONTAINER_NAME}/ && /bin/sh ${CONTAINER_NAME}.sh ${CONTAINER_NAME}

RUN rm -rf /${WORK_DIR:?} && conda clean --all --yes
