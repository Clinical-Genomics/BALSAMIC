FROM python:3.10-slim

LABEL base.image="python:3.10-slim" \
    maintainer="Clinical Genomics" \
    about.contact="support@clinicalgenomics.se" \
    software="D4tools" \
    software.version="0.3.10" \
    about.summary="Efficient storage and query format for genomic depth data" \
    about.home="https://github.com/38/d4-format" \
    about.documentation="https://github.com/38/d4-format" \
    about.license="MIT License"

ENV DEBIAN_FRONTEND=noninteractive
ENV VENV /opt/venv
ENV PATH="${VENV}/bin:$PATH"

# Install required packages and Rust toolchain
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    curl build-essential git ca-certificates liblzma-dev zlib1g-dev pkg-config \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install Rust using rustup
RUN curl https://sh.rustup.rs -sSf | bash -s -- -y
ENV PATH="/root/.cargo/bin:$PATH"

# Clone and build d4tools v0.3.10
RUN git clone --depth 1 --branch v0.3.10 https://github.com/38/d4-format.git /tmp/d4-format && \
    cd /tmp/d4-format && \
    cargo build --release && \
    cp target/release/d4tools /usr/local/bin/d4tools && \
    rm -rf /tmp/d4-format

# Optional Python venv
RUN python -m venv ${VENV} && \
    pip install --upgrade pip

# Create user and working directory
RUN adduser --disabled-password --gecos "" ubuntu && \
    chsh -s /bin/bash ubuntu && \
    mkdir -p /home/ubuntu

USER ubuntu
WORKDIR /home/ubuntu

# Default to interactive shell
CMD ["/bin/bash"]