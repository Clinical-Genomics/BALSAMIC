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

# Copy all project files to WORKDIR
COPY . .

# Execute the container specific script
WORKDIR /project/BALSAMIC/containers/${CONTAINER_NAME}
RUN /bin/sh ${CONTAINER_NAME}.sh ${CONTAINER_NAME}

# Clean up work environment
WORKDIR /
RUN rm -rf /project && conda clean --all --yes
