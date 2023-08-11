# use the ubuntu base image
FROM ubuntu:18.04

# install required packages
RUN apt-get update && apt-get install -y \
    autoconf \
    bash \
    build-essential \
    cmake \
    g++ \
    gfortran \
    git \
    libcurl4-gnutls-dev \
    hdf5-tools \
    libboost-date-time-dev \
    libboost-program-options-dev \
    libboost-system-dev \
    libboost-filesystem-dev \
    libboost-iostreams-dev \
    libbz2-dev \
    libhdf5-dev \
    libncurses-dev \
    liblzma-dev \
    zlib1g-dev \
    bcftools \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# set environment
ENV BOOST_ROOT /usr
ENV PATH="/opt/delly/bin:${PATH}"
ENV OMP_NUM_THREADS 2

# install delly
RUN cd /opt \
    && git clone --recursive https://github.com/dellytools/delly.git \
    && cd /opt/delly/ \
    && git checkout v1.0.3 \
    && make STATIC=1 PARALLEL=1 all \
    && make install \

# I'm making a change that should trigger develop containers push to Docker hub :)
