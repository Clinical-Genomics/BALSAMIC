FROM continuumio/miniconda3:4.10.3-alpine

LABEL base.image="continuumio/miniconda3:4.10.3-alpine" \
      about.home="https://github.com/Clinical-Genomics/BALSAMIC" \
      about.documentation="https://balsamic.readthedocs.io/" \
      about.license="MIT License (MIT)" \
      about.description="Bioinformatic analysis pipeline for somatic mutations in cancer"

ENV PATH="/opt/conda/bin/:${PATH}"

# Install necessary packages
RUN apk add --no-cache bash gcc git

ARG CONTAINER_NAME

# Copy the container folder into the image
COPY BALSAMIC/containers/${CONTAINER_NAME}/ /project/BALSAMIC/containers/${CONTAINER_NAME}/

# Set the working directory to the container directory
WORKDIR /project/BALSAMIC/containers/${CONTAINER_NAME}

# Run the package installation script
RUN sh ${CONTAINER_NAME}.sh ${CONTAINER_NAME}

# Install TIDDIT
WORKDIR /project
RUN python -m pip install --upgrade pip requests && \
    git clone --depth 1 --branch TIDDIT-3.8.0 https://github.com/SciLifeLab/TIDDIT.git && \
    cd TIDDIT && \
    pip install .

# Clean up work environment
WORKDIR /
RUN rm -rf /project && conda clean --all --yes
