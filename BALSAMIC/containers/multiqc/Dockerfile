FROM python:3.11.3-slim

LABEL base.image="python:3.11.3-slim"
LABEL about.home="https://github.com/Clinical-Genomics/BALSAMIC"
LABEL about.documentation="https://balsamic.readthedocs.io/"
LABEL about.license="MIT License (MIT)"
LABEL about.description="Bioinformatic analysis pipeline for somatic mutations in cancer"

# Set environment variables
ENV PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    USER=ubuntu \
    HOME=/home/ubuntu

# Install pip and multiqc, create user and home directory
RUN python -m pip install --upgrade pip && \
    pip install multiqc==1.22.3 && \
    adduser --disabled-password --gecos '' "$USER" && \
    mkdir -p "$HOME" && \
    chown -R "$USER:$USER" "$HOME"

# Switch to the non-root user
USER "$USER"
WORKDIR "$HOME"

# Set the default shell
SHELL ["/bin/bash", "-c"]

# Define the default command
CMD ["/bin/bash"]
