FROM python:3.11.3-slim

<<<<<<< HEAD
LABEL base.image="python:3.11.3-slim"
LABEL about.home="https://github.com/Clinical-Genomics/BALSAMIC"
LABEL about.documentation="https://balsamic.readthedocs.io/"
LABEL about.license="MIT License (MIT)"
LABEL about.description="Bioinformatic analysis pipeline for somatic mutations in cancer"

RUN python -m pip install --upgrade --no-cache-dir pip && \
    pip install multiqc==1.18

RUN adduser --disabled-password --gecos '' ubuntu && \
chsh -s /bin/bash && \
mkdir -p /home/ubuntu

USER ubuntu
WORKDIR /home/ubuntu
=======
LABEL base.image="python:3.11.3-slim" \
    maintainer="Clinical Genomics" \
    about.contact="support@clinicalgenomics.se" \
    about.home="https://github.com/MultiQC/MultiQC" \
    software.version="1.22.3" \
    about.documentation="https://multiqc.info/" \
    about.license="MIT License (MIT)" \
    about.description="Aggregate bioinformatics results across many samples into a single report"

# Set environment variables
ENV PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    USER=ubuntu \
    HOME=/home/ubuntu

WORKDIR "$HOME"

# Install pip and multiqc, create user and home directory
RUN python -m pip install --upgrade pip && \
    pip install multiqc==1.22.3 && \
    adduser --disabled-password --gecos '' "$USER" && \
    chown -R "$USER:$USER" "$HOME"

# Switch to the non-root user
USER "$USER"

# Set the default shell
SHELL ["/bin/bash", "-c"]

# Define the default command
>>>>>>> develop
CMD ["/bin/bash"]
