FROM markusriester/purecn:2.6.4-amd64

LABEL base.image="markusriester/purecn:2.6.4-amd64"
LABEL maintainer="Clinical Genomics"
LABEL about.contact="support@clinicalgenomics.se"
LABEL software="PureCN"
LABEL software.version="2.6.4"
LABEL about.summary="Copy number calling and variant classification using targeted short read sequencing"
LABEL about.home="https://github.com/lima1/PureCN"
LABEL about.documentation="https://bioconductor.org/packages/devel/bioc/html/PureCN.html"
LABEL about.license="MIT License (MIT)"

RUN apt-get update && apt-get -y upgrade && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN adduser --disabled-password --gecos '' ubuntu && \
    chsh -s /bin/bash && mkdir -p /home/ubuntu

USER ubuntu
WORKDIR /home/ubuntu
CMD ["/bin/bash"]
