FROM alpine:3.11.5

LABEL base.image="alpine:3.11.5" \
      maintainer="Clinical Genomics" \
      about.contact="support@clinicalgenomics.se" \
      software="somalier" \
      software.version="0.2.19" \
      about.summary="Fast sample-swap and relatedness checks on BAMs/CRAMs/VCFs/GVCFs" \
      about.home="https://github.com/brentp/somalier" \
      about.documentation="https://github.com/brentp/somalier/wiki" \
      about.license="MIT License (MIT)"

ARG nim_version=1.6.20
ENV CFLAGS="-fPIC -O3"
ENV HTTPS_METHOD=noredirect

RUN apk upgrade && apk add --no-cache wget git xz bzip2-static musl m4 \
    autoconf tar xz-dev bzip2-dev build-base libpthread-stubs libzip-dev gfortran \
    openssl-libs-static openblas-static pcre-dev curl llvm-dev curl-static \
    bash curl-dev clang-static nghttp2-static  zlib-static cmake procps

RUN mkdir -p /usr/local/include && \
    git clone -b v1.18 --depth 1 https://github.com/ebiggers/libdeflate.git && \
    cd libdeflate && cmake -B build && cmake --build build && cmake --install build  && \
    ln -s /usr/local/lib64/* /usr/local/lib && \
    cd .. && rm -rf libdeflate

RUN cd / && \
    wget --max-redirect=0 -q https://nim-lang.org/download/nim-${nim_version}-linux_x64.tar.xz && \
    tar xf nim-${nim_version}-linux_x64.tar.xz && \
    echo 'PATH=/nim-${nim_version}/bin:$PATH' >> ~/.bashrc && \
    echo 'PATH=/nim-${nim_version}/bin:$PATH' >> ~/.bash_profile && \
    echo 'PATH=/nim-${nim_version}/bin:$PATH' >> /etc/environment  && \
    rm -f nim-${nim_version}-linux_x64.tar.xz

RUN apk add --no-cache cmake openssl-dev && \
    wget --max-redirect=0 -q https://libzip.org/download/libzip-1.6.1.tar.gz && \
	tar xzvf libzip-1.6.1.tar.gz && \
	cd libzip-1.6.1 && \
	mkdir build && cd build && \
	cmake -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=/usr/local/ ../ && \
	make -j4 CFLAGS="-fPIC -O3" install && \
	cd ../../ && rm -rf libzip-1.6.1*

ENV PATH=:/root/.nimble/bin:/nim-${nim_version}/bin/:$PATH

RUN git clone --depth 1 -b 1.18 --recursive https://github.com/samtools/htslib && \
    cd htslib && autoheader && autoconf && \
    ./configure --enable-s3 --enable-gcs --enable-libcurl --with-libdeflate && \
    make -j4 CFLAGS="-fPIC -O3" install && \
    cd ../ && \
    git clone --depth 1 -b 1.18 --recursive https://github.com/samtools/bcftools && \
    cd bcftools && autoheader && autoconf && \
    ./configure --enable-s3 --enable-libcurl --with-libdeflate && \
    make -j4 CFLAGS="-fPIC -O3" install && \
    cd ../ && rm -rf htslib bcftools

RUN cd / && git clone -b v0.2.19 --depth 5 https://github.com/brentp/somalier.git && \
    nimble install -y nimble@0.14.2 && \
    cd somalier && \
    sed -i 's/"arraymancer"/"arraymancer <= 0.7.32"/' somalier.nimble && \
    nimble install -d -y

RUN apk add --no-cache openblas && \
    cd /somalier && \
    nim c -d:danger -d:nsb_static -d:release -d:openmp -d:blas=openblas -d:lapack=openblas -o:/usr/bin/somalier src/somalier && \
    cp scripts/ancestry-labels-1kg.tsv / && \
    rm -rf /somalier && somalier --help

ENV somalier_ancestry_labels=/ancestry_labels-1kg.tsv

RUN rm -rf /var/cache/apk/*

RUN adduser -D defaultuser

USER defaultuser
