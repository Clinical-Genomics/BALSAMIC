FROM ubuntu:22.04

LABEL base.image="ubuntu:22.04" \
      maintainer="Clinical Genomics" \
      about.contact="support@clinicalgenomics.se" \
      software="msisensor-pro" \
      software.version="1.2.0" \
      about.summary="Microsatellite Instability (MSI) detection using high-throughput sequencing data" \
      about.home="https://github.com/xjtu-omics/msisensor-pro" \
      about.documentation="https://github.com/xjtu-omics/msisensor-pro/wiki" \
      about.license="MIT License (MIT)"

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && apt-get -y upgrade && \
    apt-get -y install -y --no-install-recommends git make build-essential \
    liblzma-dev libbz2-dev zlib1g-dev libncurses5-dev libncursesw5-dev ca-certificates && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN cd /tmp && \
  git clone https://github.com/xjtu-omics/msisensor-pro && \
  cd msisensor-pro && \
  /tmp/msisensor-pro/INSTALL && \
  cp -r /tmp/msisensor-pro/binary/msisensor-pro /usr/bin/ && \
  cp -r /tmp/msisensor-pro/vendor/htslib-1.11/libhts.so* /lib/x86_64-linux-gnu && \
  rm -rf /tmp/msisensor-pro

RUN adduser --disabled-password --gecos '' ubuntu && \
    chsh -s /bin/bash && mkdir -p /home/ubuntu

USER ubuntu
WORKDIR /home/ubuntu
CMD ["/bin/bash"]
