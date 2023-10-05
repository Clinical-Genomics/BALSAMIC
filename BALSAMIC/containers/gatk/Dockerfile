FROM python:3.11-slim

LABEL base.image="python:3.11-slim"
LABEL maintainer="Clinical Genomics"
LABEL about.contact="support@clinicalgenomics.se"
LABEL software="GATK"
LABEL software.version="4.4.0.4"
LABEL about.summary="Variant Discovery in High-Throughput Sequencing Data"
LABEL about.home="https://github.com/broadinstitute/gatk"
LABEL about.documentation="https://gatk.broadinstitute.org/hc/en-us"
LABEL about.license="MIT License (MIT)"

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && apt-get -y upgrade && \
    apt-get -y install --no-install-recommends wget openjdk-17-jre-headless apt-utils \
    git r-base-core && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV JAVA_HOME=/usr/lib/jvm/java-17-openjdk-amd64

WORKDIR /tmp
# install GATK4
ENV GATK_VERSION="4.4.0.0"
RUN wget --no-verbose https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip && \
    unzip gatk-${GATK_VERSION}.zip -d /opt && \
    rm gatk-${GATK_VERSION}.zip

ENV PATH /opt/gatk-${GATK_VERSION}:$PATH

RUN adduser --disabled-password --gecos '' ubuntu && \
    chsh -s /bin/bash && mkdir -p /home/ubuntu

USER ubuntu
WORKDIR /home/ubuntu
CMD ["/bin/bash"]

