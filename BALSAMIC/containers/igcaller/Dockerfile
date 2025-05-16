# Base image with Python 3
FROM python:3.8-slim

LABEL maintainer="Your Name <your.email@example.com>"
LABEL description="Docker image for IgCaller v1.3 with samtools v1.9"

ENV DEBIAN_FRONTEND=noninteractive \
    IG_CALLER_DIR=/opt/IgCaller \
    SAMTOOLS_VERSION=1.9

# Install dependencies for building samtools
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    make \
    gcc \
    libbz2-dev \
    liblzma-dev \
    libncurses-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    wget \
    unzip \
    ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install samtools v1.9
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && \
    ./configure --prefix=/usr/local && \
    make && make install && \
    cd .. && rm -rf samtools-${SAMTOOLS_VERSION} samtools-${SAMTOOLS_VERSION}.tar.bz2

# Clone IgCaller at release v1.3 and set permissions
RUN git clone --branch v1.3 https://github.com/ferrannadeu/IgCaller.git $IG_CALLER_DIR && \
    chmod +x $IG_CALLER_DIR/IgCaller

# Set working directory
WORKDIR $IG_CALLER_DIR

# Install Python dependencies (versions tested by authors)
RUN pip install --no-cache-dir \
    regex==2020.11.13 \
    argparse==1.4.0 \
    numpy==1.16.3 \
    scipy==1.3.0

# Add IgCaller to PATH
ENV PATH="$IG_CALLER_DIR:$PATH"
