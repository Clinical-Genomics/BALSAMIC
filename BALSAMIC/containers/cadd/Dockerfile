FROM condaforge/mambaforge:24.9.2-0

LABEL base.image="condaforge/miniconda3:25.1.1-2"
LABEL maintainer="Clinical Genomics"
LABEL about.contact="support@clinicalgenomics.se"
LABEL software="cadd"
LABEL software.version="1.6"
LABEL about.summary="Combined Annotation Dependent Depletion (CADD)"
LABEL about.home="https://cadd.gs.washington.edu/"
LABEL about.documentation="https://cadd.gs.washington.edu/info"
LABEL about.license="MIT License (MIT)"

## Install snakemake and unzip. Download and extract CADD repository
RUN conda install -y -c conda-forge -c bioconda snakemake=7.32.2 unzip && \
    conda clean -afy && \
    wget --no-verbose \
    https://github.com/kircherlab/CADD-scripts/archive/CADD1.6.zip \
    -O /opt/conda/share/CADD1.6.zip && \
    unzip -o /opt/conda/share/CADD1.6.zip -d /opt/conda/share/ && \
    rm /opt/conda/share/CADD1.6.zip && \
    mv /opt/conda/share/CADD-scripts-CADD1.6 /opt/conda/share/CADD-scripts

## Install CADD-scripts
RUN cd /opt/conda/share/CADD-scripts && \
    snakemake test/input.tsv.gz --use-conda --conda-create-envs-only \
    --conda-prefix /opt/conda/share/CADD-scripts/envs \
    --cores 1 --configfile config/config_GRCh38_v1.6.yml --snakefile Snakefile && \
    ln -s /opt/conda/share/CADD-scripts/CADD.sh /opt/conda/bin/CADD.sh && \
    conda clean -yfa

# Run as a non-root
RUN useradd -m dockeruser
USER dockeruser

# Set the work directory
WORKDIR /home/dockeruser/

# Set the default shell
SHELL ["/bin/bash", "-c"]

# Define the default command
CMD ["/bin/bash"]
