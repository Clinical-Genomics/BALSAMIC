FROM continuumio/miniconda3:25.3.1-1

LABEL base.image="continuumio/miniconda3" \
    maintainer="Clinical Genomics" \
    about.contact="support@clinicalgenomics.se" \
    about.license="MIT License (MIT)" \
    about.home="https://github.com/Clinical-Genomics/BALSAMIC" \
    about.documentation="https://balsamic.readthedocs.io/" \
    about.description="Bioinformatic analysis pipeline for somatic mutations in cancer"

ARG CONTAINER_NAME

# Set environment variables
ENV GATK_VERSION=4.6.2

# Create a conda environment for GATK4
# Disable conda auto-updates, which can cause compatibility issues with GATK
# Enable the (much) faster libmamba solver
RUN conda config --set auto_update_conda false && \
    conda config --set solver libmamba && \
    conda install -y -c bioconda -c conda-forge gatk4=${GATK_VERSION} && \
    conda clean -afy

#ENV PATH="/opt/conda/envs/${CONTAINER_NAME}/bin/:${PATH}"
ENV PATH="/opt/conda/bin:${PATH}"

# Run as a non-root
RUN useradd -m dockeruser
USER dockeruser

# Set the work directory
WORKDIR /home/dockeruser/

# Set the default shell
SHELL ["/bin/bash", "-c"]

# Define the default command
CMD ["/bin/bash"]
