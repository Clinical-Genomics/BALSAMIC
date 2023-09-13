FROM python:3.10-slim


LABEL base.image="python:3.10-slim"
LABEL maintainer="Clinical Genomics"
LABEL about.contact="support@clinicalgenomics.se"
LABEL software="cnvkit"
LABEL software.version="0.9.10"
LABEL about.summary="Copy number variant detection from targeted DNA sequencing"
LABEL about.home="https://github.com/etal/cnvkit"
LABEL about.documentation="https://cnvkit.readthedocs.io"
LABEL about.license="MIT License (MIT)"

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && apt-get -y upgrade && \
    apt-get -y install --no-install-recommends tabix liblzma-dev zlib1g-dev \
    r-base-core r-bioc-dnacopy && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV VENV /opt/venv

RUN python -m pip install --upgrade --no-cache-dir pip
RUN python -m venv ${VENV}
ENV PATH="${VENV}/bin:$PATH"

RUN pip install --no-cache-dir --upgrade pip

RUN pip install --no-cache-dir cnvkit==0.9.10

RUN adduser --disabled-password --gecos '' ubuntu && \
    chsh -s /bin/bash && mkdir -p /home/ubuntu

USER ubuntu
WORKDIR /home/ubuntu
CMD ["/bin/bash"]
