FROM continuumio/miniconda3:25.1.1-2

LABEL base.image="continuumio/miniconda3" \
    software="cadd" \
    software.version="1.6" \
    about.summary="Combined Annotation Dependent Depletion (CADD)" \
    about.home="https://cadd.gs.washington.edu/" \
    about.documentation="https://cadd.gs.washington.edu/info" \
    maintainer="Clinical Genomics" \
    about.contact="support@clinicalgenomics.se" \
    about.license="MIT License (MIT)"

ARG CONTAINER_NAME

# Create, download and create conda environments
RUN conda config --set channel_priority strict && \
    conda install -c conda-forge -c bioconda 'snakemake=8' unzip && \
    conda clean -afy  && \
    wget --no-verbose \
    https://github.com/kircherlab/CADD-scripts/archive/refs/tags/v1.7.2.zip \
    -O /opt/conda/share/CADD.zip && \
    unzip -o /opt/conda/share/CADD.zip -d /opt/conda/share/ && \
    rm /opt/conda/share/CADD.zip && \
    export CADD=/opt/conda/share/CADD-scripts-1.7.2 && \
    cd $CADD  && \
    snakemake test/input.tsv.gz \
    --software-deployment-method conda \
    --conda-create-envs-only \
    --conda-prefix envs/conda \
    --configfile config/config_GRCh38_v1.7.yml \
    --snakefile Snakefile -c 1 && \
    ln -s $CADD/CADD.sh /opt/conda/bin/CADD.sh

# Set CADD as enviromental variable
ENV CADD=/opt/conda/share/CADD-scripts-1.7.2

# Run as a non-root
RUN useradd -m dockeruser
USER dockeruser

# Set the work directory
WORKDIR /home/dockeruser/

# Set the default shell
SHELL ["/bin/bash", "-c"]

# Define the default command
CMD ["/bin/bash"]
