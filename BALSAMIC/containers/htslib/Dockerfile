FROM ubuntu:22.04

LABEL base.image="ubuntu:22.04"
LABEL maintainer="Clinical Genomics"
LABEL about.contact="support@clinicalgenomics.se"
LABEL software="htslib, samtools, bcftools, tabix"
LABEL software.version="1.13"
LABEL about.summary="A unified C library for accessing common bioinformatics file formats"
LABEL about.home="http://www.htslib.org/"
LABEL about.documentation="http://www.htslib.org/doc/#manual-pages"
LABEL about.license="MIT License (MIT)"

RUN apt-get update && \
     apt-get -y upgrade && \
     apt-get -y  install --no-install-recommends \
     samtools=1.13-4 bcftools=1.13-1 tabix=1.13+ds-2build1 && \
     apt-get clean && \
     rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN adduser --disabled-password --gecos '' ubuntu &&\
chsh -s /bin/bash && mkdir -p /home/ubuntu

USER    ubuntu
WORKDIR /home/ubuntu
CMD ["/bin/bash"]
