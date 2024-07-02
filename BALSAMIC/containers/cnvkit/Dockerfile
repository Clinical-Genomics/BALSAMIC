FROM python:3.10-slim

LABEL base.image="python:3.10-slim" \
    maintainer="Clinical Genomics" \
    about.contact="support@clinicalgenomics.se" \
    software="CNVkit" \
    software.version="0.9.10" \
    about.summary="Copy number variant detection from targeted DNA sequencing" \
    about.home="https://github.com/etal/cnvkit" \
    about.documentation="https://cnvkit.readthedocs.io" \
    about.license="MIT License (MIT)"

ENV DEBIAN_FRONTEND noninteractive
ENV VENV /opt/venv
ENV PATH="${VENV}/bin:$PATH"

RUN apt-get update && apt-get -y upgrade && \
    apt-get -y install --no-install-recommends tabix liblzma-dev zlib1g-dev \
    r-base-core r-bioc-dnacopy && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN python -m pip install --upgrade --no-cache-dir pip && \
    python -m venv ${VENV} && \
    pip install --no-cache-dir "cnvkit==0.9.10" "numpy<2.0.0"

RUN adduser --disabled-password --gecos "" ubuntu && \
    chsh -s /bin/bash ubuntu && \
    mkdir -p /home/ubuntu

USER ubuntu
WORKDIR /home/ubuntu
CMD ["/bin/bash"]
