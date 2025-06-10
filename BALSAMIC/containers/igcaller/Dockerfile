FROM python:3.8-slim

LABEL maintainer="blabla"
LABEL description="IgCaller v1.3 with Samtools v1.9"

ENV DEBIAN_FRONTEND=noninteractive
ENV IG_CALLER_DIR=/opt/IgCaller
ENV SAMTOOLS_VERSION=1.9
ENV HTSLIB_VERSION=1.9

# Install build tools and dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    gfortran \
    autoconf \
    automake \
    libtool \
    pkg-config \
    libatlas-base-dev \
    libffi-dev \
    libbz2-dev \
    liblzma-dev \
    libssl-dev \
    libsqlite3-dev \
    libreadline-dev \
    libncursesw5-dev \
    libxml2-dev \
    libxmlsec1-dev \
    zlib1g-dev \
    libblas-dev \
    liblapack-dev \
    wget \
    git \
    ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# Download and install htslib (required for samtools)
RUN wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 && \
    tar -xjf htslib-${HTSLIB_VERSION}.tar.bz2 && \
    cd htslib-${HTSLIB_VERSION} && \
    ./configure --prefix=/usr/local && make && make install && \
    cd .. && rm -rf htslib-${HTSLIB_VERSION}*

# Download and install samtools v1.9
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && \
    ./configure --prefix=/usr/local && make && make install && \
    cd .. && rm -rf samtools-${SAMTOOLS_VERSION}*

# Clone IgCaller at branch v1.3 and set permissions
RUN git clone https://github.com/ferrannadeu/IgCaller.git /opt/IgCaller && \
    cd /opt/IgCaller && \
    git fetch --tags && \
    git checkout tags/v1.3

# Set working directory
WORKDIR $IG_CALLER_DIR

# Install wheel to enable binary wheels
RUN pip install --upgrade pip setuptools wheel

# Install compatible versions separately
RUN pip install --no-cache-dir \
    numpy==1.18.5 \
    regex==2020.11.13

# Then install scipy (requires numpy to be pre-installed)
RUN pip install --no-cache-dir \
    scipy==1.4.1

# Create a wrapper script that runs IgCaller
RUN echo '#!/bin/bash\npython3 /opt/IgCaller/IgCaller.py "$@"' > /usr/local/bin/igcaller && \
    chmod +x /usr/local/bin/igcaller