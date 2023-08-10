FROM condaforge/mambaforge:23.1.0-4

LABEL base.image="condaforge/mambaforge:23.1.0.4"
LABEL about.home="https://github.com/Clinical-Genomics/BALSAMIC"
LABEL about.documentation="https://balsamic.readthedocs.io/"
LABEL about.license="MIT License (MIT)"
LABEL about.description="CADD for BALSAMIC"

RUN apt-get update && \
     apt-get -y upgrade && \
     apt-get -y install -y --no-install-recommends \
     tabix \
     unzip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

## Install snakemake
RUN conda install -c conda-forge -c bioconda snakemake=7.32.2

## Donwload CADD repository
RUN wget --no-verbose \
    https://github.com/kircherlab/CADD-scripts/archive/CADD1.6.zip \
    -O /opt/conda/share/CADD1.6.zip

## Extract CADD directory
RUN unzip -o /opt/conda/share/CADD1.6.zip -d /opt/conda/share/ && \
    rm /opt/conda/share/CADD1.6.zip && \
    mv /opt/conda/share/CADD-scripts-CADD1.6 /opt/conda/share/CADD-scripts

## Install CADD-scripts
RUN cd /opt/conda/share/CADD-scripts && \
    snakemake test/input.tsv.gz --use-conda --conda-create-envs-only \
    --conda-prefix /opt/conda/share/CADD-scripts/envs \
    --cores 1 --configfile config/config_GRCh38_v1.6.yml --snakefile Snakefile
RUN  ln -s /opt/conda/share/CADD-scripts/CADD.sh /opt/conda/bin/CADD.sh

## Clean-up Conda
RUN /opt/conda/bin/conda clean -ya

RUN adduser --disabled-password --gecos '' ubuntu &&\
chsh -s /bin/bash && mkdir -p /home/ubuntu

USER    ubuntu
WORKDIR /home/ubuntu
CMD ["/bin/bash"]
