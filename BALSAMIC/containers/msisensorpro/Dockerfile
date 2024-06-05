FROM ubuntu:22.04

LABEL base.image="ubuntu:22.04"
LABEL maintainer="Clinical Genomics"
LABEL about.contact="support@clinicalgenomics.se"
LABEL software="msisensor-pro"
LABEL software.version="1.2.0"
LABEL about.summary="Microsatellite Instability (MSI) detection using high-throughput sequencing data"
LABEL about.home="https://github.com/xjtu-omics/msisensor-pro"
LABEL about.documentation="https://github.com/xjtu-omics/msisensor-pro/wiki"
LABEL about.license="MIT License (MIT)"

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && apt-get -y upgrade && \
    apt-get -y install -y  --no-install-recommends git make build-essential \
    liblzma-dev libbz2-dev zlib1g-dev libncurses5-dev libncursesw5-dev ca-certificates && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN cd /tmp && \
  git clone https://github.com/xjtu-omics/msisensor-pro && \
  cd msisensor-pro && \
  git checkout v1.2.0 && \
  ./INSTALL

RUN  cp -r /tmp/msisensor-pro/binary/msisensor-pro /usr/bin/ && \
  cp -r /tmp/msisensor-pro/vendor/htslib-1.11/libhts.so* /lib/x86_64-linux-gnu

RUN adduser --disabled-password --gecos '' ubuntu && \
    chsh -s /bin/bash && mkdir -p /home/ubuntu

USER ubuntu
WORKDIR /home/ubuntu
CMD ["/bin/bash"]
