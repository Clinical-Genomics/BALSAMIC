FROM python:3.11.3-slim AS builder

ARG CONTAINER_NAME
ENV WORK_DIR /opt/
ENV VENV ${WORK_DIR}/venv

RUN apt-get update && \
    apt-get -y upgrade && \
    apt-get -y install --no-install-recommends git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN python -m pip install --upgrade --no-cache-dir pip && \
    python -m venv ${VENV}

ENV PATH="${VENV}/bin:$PATH"

RUN cd ${WORK_DIR} && \
    git clone https://github.com/NBISweden/vcf2cytosure.git && \
    cd ${WORK_DIR}/${CONTAINER_NAME} && \
    git checkout v0.9.1 && \
    pip install --no-cache-dir .

FROM python:3.11-slim

LABEL base.image="python:3.11-slim" \
      about.home="https://github.com/Clinical-Genomics/BALSAMIC" \
      about.documentation="https://balsamic.readthedocs.io/" \
      about.license="MIT License (MIT)" \
      about.description="vcf2cytosure for BALSAMIC"

RUN apt-get update && \
    apt-get -y upgrade && \
    apt-get -y install --no-install-recommends tabix && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV WORK_DIR /opt/
ENV VENV ${WORK_DIR}/venv
ENV PATH="${VENV}/bin:$PATH"

COPY --from=builder ${WORK_DIR} ${WORK_DIR}

RUN adduser --disabled-password --gecos '' ubuntu && \
    chsh -s /bin/bash ubuntu && \
    mkdir -p /home/ubuntu

USER ubuntu
WORKDIR /home/ubuntu
CMD ["/bin/bash"]
