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

# Set the working directory
WORKDIR /project

# Copy the relevant container folder into the image
COPY BALSAMIC/containers/${CONTAINER_NAME}/ /project/BALSAMIC/containers/${CONTAINER_NAME}/

# Set the working directory to the container directory
WORKDIR /project/BALSAMIC/containers/${CONTAINER_NAME}

# Run the script
RUN sh ${CONTAINER_NAME}.sh ${CONTAINER_NAME}

# Clean up work environment
WORKDIR /
RUN rm -rf /project && conda clean --all --yes
