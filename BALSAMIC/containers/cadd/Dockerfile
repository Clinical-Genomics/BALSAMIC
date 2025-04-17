FROM continuumio/miniconda3

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
COPY BALSAMIC/containers/${CONTAINER_NAME}/${CONTAINER_NAME}.yaml environment.yaml

# Create conda environment
RUN conda env create -n ${CONTAINER_NAME} -f environment.yaml && \
    conda clean -afy && \
    rm -f environment.yaml

ENV PATH="/opt/conda/envs/${CONTAINER_NAME}/bin/:${PATH}"

# The cadd-scripts conda package is available as cadd.sh instead
# of CADD.sh. Make link to keep compatibility
RUN ln -s /opt/conda/envs/cadd/bin/cadd.sh /opt/conda/envs/cadd/bin/CADD.sh

# Set the default shell
SHELL ["/bin/bash", "-c"]

# Define the default command
CMD ["/bin/bash"]
