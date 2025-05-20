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

# Copy enviroment file
COPY BALSAMIC/containers/${CONTAINER_NAME}/${CONTAINER_NAME}.yaml /tmp/environment.yaml

# Create conda environment
# The cadd-scripts conda package is available as cadd.sh instead
# of CADD.sh so we make link to keep compatibility
RUN conda env create -n ${CONTAINER_NAME} -f /tmp/environment.yaml && \
    conda clean -afy && \
    rm -f /tmp/environment.yaml && \
    ln -s /opt/conda/envs/cadd/bin/cadd.sh /opt/conda/envs/cadd/bin/CADD.sh && \
    mkdir -p /opt/conda/envs/cadd/share/cadd-scripts-1.7.2-1/envs/apptainer && \
    mkdir -p /opt/conda/envs/cadd/share/cadd-scripts-1.7.2-1/envs/conda

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
