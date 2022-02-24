FROM continuumio/miniconda3:4.9.2-alpine

LABEL base_image="continuumio/miniconda3:4.9.2-alpine"
LABEL about.home="https://github.com/Clinical-Genomics/BALSAMIC"
LABEL about.documentation="https://balsamic.readthedocs.io/"
LABEL about.license="MIT License (MIT)"
LABEL about.maintainer="Ashwini Jeggari ashwini dot jeggari at scilifelab dot se"
LABEL about.description="Bioinformatic analysis pipeline for somatic mutations in cancer"

ENV PATH="/opt/vcf2cytosure/bin:${PATH}"

RUN apk add --no-cache bash gcc git python3

RUN cd /opt \
    && git clone https://github.com/NBISweden/vcf2cytosure.git \
    && cd vcf2cytosure/ \
    && pip install -e .
