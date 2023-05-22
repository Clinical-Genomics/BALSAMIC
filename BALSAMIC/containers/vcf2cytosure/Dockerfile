FROM python:3.11-slim AS builder

ARG CONTAINER_NAME
ENV WORK_DIR /opt/
ENV VENV ${WORK_DIR}/venv

RUN apt-get update && \
    apt-get -y upgrade && \
    apt-get -y install -y --no-install-recommends git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN python -m pip install --upgrade --no-cache-dir pip
RUN python -m venv ${VENV}
ENV PATH="${VENV}/bin:$PATH"

RUN pip install --no-cache-dir \
cyvcf2==0.30.22 \
lxml==4.9.2 \
pandas==2.0.1

RUN cd /opt &&\
    git clone https://github.com/NBISweden/vcf2cytosure.git && \
    cd /opt/${CONTAINER_NAME}/ && \
    git checkout v0.8 && \
	pip install --no-cache-dir -e .

from python:3.11-slim

LABEL base.image="python:3.11-slim"
LABEL about.home="https://github.com/Clinical-Genomics/BALSAMIC"
LABEL about.documentation="https://balsamic.readthedocs.io/"
LABEL about.license="MIT License (MIT)"
LABEL about.description="vcf2cytosure for BALSAMIC"

RUN apt-get update && \
    apt-get -y upgrade && \
    apt-get -y install -y --no-install-recommends tabix && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV WORK_DIR /opt/
ENV VENV ${WORK_DIR}/venv

COPY --from=builder ${WORK_DIR} ${WORK_DIR}

ENV PATH="${VENV}/bin:$PATH"

RUN adduser --disabled-password --gecos '' ubuntu && \
chsh -s /bin/bash && \
mkdir -p /home/ubuntu

USER ubuntu
WORKDIR /home/ubuntu
CMD ["/bin/bash"]
