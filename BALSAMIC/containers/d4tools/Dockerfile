# Use a lightweight base image with support for conda
FROM continuumio/miniconda3:latest

# Install system dependencies
RUN apt-get update && \
    apt-get install -y git build-essential && \
    rm -rf /var/lib/apt/lists/*

# Create conda environment and install a specific d4tools version
RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda install -y d4tools=0.3.10 && \
    conda clean --all -y

# Verify installation
RUN d4tools --version

# Set entrypoint (optional)
ENTRYPOINT ["d4tools"]