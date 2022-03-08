FROM continuumio/miniconda3:4.10.3-alpine

LABEL base.image="continuumio/miniconda3:4.10.3-alpine"
LABEL about.home="https://github.com/Clinical-Genomics/BALSAMIC"
LABEL about.documentation="https://balsamic.readthedocs.io/"
LABEL about.license="MIT License (MIT)"
LABEL about.maintainer="Ashwini Jeggari ashwini dot jeggari at scilifelab dot se"
LABEL about.description="Bioinformatic analysis pipeline for somatic mutations in cancer"

ARG CONTAINER_NAME
ENV PATH="/opt/${CONTAINER_NAME}/bin:${PATH}"
ENV PYTHONPATH="/opt/${CONTAINER_NAME}"


RUN apk add --no-cache bash gcc git python3

RUN cd /opt \
    && git clone https://github.com/NBISweden/vcf2cytosure.git \
    && cd /opt/${CONTAINER_NAME}/  \
    && pip install --no-cache-dir .
