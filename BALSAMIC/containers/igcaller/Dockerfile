FROM python:3.8-slim

LABEL maintainer="Your Name <your.email@example.com>"
LABEL description="IgCaller v1.3 with Samtools v1.9"

ENV DEBIAN_FRONTEND=noninteractive
ENV IG_CALLER_DIR=/opt/IgCaller
ENV SAMTOOLS_VERSION=1.9
ENV HTSLIB_VERSION=1.9

# Install build tools and dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    build-essential \
    wget \
    curl \
    make \
    gcc \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    zlib1g-dev \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

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
    git checkout -b igcaller-v1.3

# Set working directory
WORKDIR $IG_CALLER_DIR

# Install Python dependencies
RUN pip install --no-cache-dir \
    regex==2020.11.13 \
    numpy==1.16.3 \
    scipy==1.3.0

# Add IgCaller to PATH
ENV PATH="${IG_CALLER_DIR}:$PATH"
