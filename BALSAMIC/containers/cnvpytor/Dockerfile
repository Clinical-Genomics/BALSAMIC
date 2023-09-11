FROM python:3.11-slim AS builder

LABEL base.image="python:3.11-slim"
LABEL maintainer="Clinical Genomics"
LABEL about.contact="support@clinicalgenomics.se"
LABEL software="cnvpytor"
LABEL software.version="1.3.1"
LABEL about.summary="A tool for CNV discovery and genotyping from depth-of-coverage by mapped reads"
LABEL about.home="https://github.com/abyzovlab/CNVpytor"
LABEL about.documentation="https://github.com/abyzovlab/CNVpytor/wiki"
LABEL about.license="MIT License (MIT)"

ARG CONTAINER_NAME
ENV WORK_DIR /opt
ENV VENV ${WORK_DIR}/venv

RUN apt-get update && \
     apt-get -y upgrade && \
     apt-get -y install -y --no-install-recommends git && \
     apt-get clean && \
     rm -rf /var/lib/apt/lists/*

RUN python -m pip install --upgrade pip
RUN python -m venv ${VENV}
ENV PATH="${VENV}/bin:$PATH"

RUN pip install seaborn==0.12.2

RUN cd /opt &&\
    git clone https://github.com/abyzovlab/CNVpytor.git &&\
    cd /opt/${CONTAINER_NAME}/ &&\
    git checkout v1.3.1 &&\
	pip install --no-cache-dir -e . && \
    cnvpytor -download

from python:3.11-slim

RUN apt-get update && \
     apt-get -y upgrade && \
     apt-get -y install -y --no-install-recommends tabix && \
     apt-get clean && \
     rm -rf /var/lib/apt/lists/*

ENV WORK_DIR /opt
ENV VENV ${WORK_DIR}/venv

COPY --from=builder /opt /opt

ENV PATH="${VENV}/bin:$PATH"

RUN adduser --disabled-password --gecos '' ubuntu &&\
chsh -s /bin/bash && mkdir -p /home/ubuntu

USER    ubuntu
WORKDIR /home/ubuntu
CMD ["/bin/bash"]
