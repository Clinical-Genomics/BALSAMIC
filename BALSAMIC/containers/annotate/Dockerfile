FROM continuumio/miniconda3:25.1.1-2

LABEL base.image="continuumio/miniconda3" \
    maintainer="Clinical Genomics" \
    about.contact="support@clinicalgenomics.se" \
    about.license="MIT License (MIT)" \
    about.home="https://github.com/Clinical-Genomics/BALSAMIC" \
    about.documentation="https://balsamic.readthedocs.io/" \
    about.description="Bioinformatic analysis pipeline for somatic mutations in cancer"

ARG CONTAINER_NAME

# Copy enviroment file
COPY BALSAMIC/containers/${CONTAINER_NAME}/${CONTAINER_NAME}.yaml /tmp/environment.yaml

# Create conda environment
RUN conda env create -n ${CONTAINER_NAME} -f /tmp/environment.yaml && \
    conda clean -afy && \
    rm -f /tmp/environment.yaml

ENV PATH="/opt/conda/envs/${CONTAINER_NAME}/bin/:${PATH}"

# Run as a non-root
RUN useradd -m dockeruser
USER dockeruser

# Set the work directory
WORKDIR /home/dockeruser/

# Set the default shell
SHELL ["/bin/bash", "-c"]

# Define the default command
CMD ["/bin/bash"]
