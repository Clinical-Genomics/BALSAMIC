FROM quay.io/wtsicgp/pcap-core:5.4.2 as builder

USER root

ENV VER_CGPVCF="v2.2.1"
ENV VER_VCFTOOLS="0.1.16"
ENV VER_ALLELECOUNT="v4.2.1"

RUN apt-get -yq update
RUN apt-get install -qy --no-install-recommends lsb-release
RUN apt-get install -qy --no-install-recommends gnupg
RUN apt-get install -yq --no-install-recommends dirmngr
RUN apt-get install -yq --no-install-recommends wget
RUN apt-get install -yq --no-install-recommends build-essential
RUN apt-get install -yq --no-install-recommends ca-certificates
RUN apt-get install -yq --no-install-recommends locales
RUN apt-get install -yq --no-install-recommends g++
RUN apt-get install -yq --no-install-recommends make
RUN apt-get install -yq --no-install-recommends gcc
RUN apt-get install -yq --no-install-recommends git
ENV DEBIAN_FRONTEND "noninteractive"
RUN apt-get -yq update
RUN apt-get install -yq --no-install-recommends software-properties-common
RUN echo "deb http://cran.rstudio.com/bin/linux/ubuntu $(lsb_release -cs)-cran40/" >> /etc/apt/sources.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
RUN apt-get -yq update
ENV DEBIAN_FRONTEND "noninteractive"
RUN apt-get install -yq --no-install-recommends locales
RUN apt-get install -yq --no-install-recommends g++
RUN apt-get install -yq --no-install-recommends make
RUN apt-get install -yq --no-install-recommends gcc
RUN apt-get install -yq --no-install-recommends wget
RUN apt-get install -yq --no-install-recommends pkg-config
RUN apt-get install -yq --no-install-recommends zlib1g-dev
RUN apt-get install -yq --no-install-recommends libbz2-dev
RUN apt-get install -yq --no-install-recommends unzip
RUN apt-get install -yq --no-install-recommends libpng-dev
RUN apt-get install -yq --no-install-recommends tzdata
RUN apt-get install -yq --no-install-recommends libcurl4-openssl-dev
RUN apt-get install -yq --no-install-recommends libxml2-dev
RUN apt-get install -yq --no-install-recommends libgit2-dev
RUN apt-get install -yq --no-install-recommends liblzma-dev
RUN apt-get install -yq --no-install-recommends libssl-dev
RUN apt-get install -yq --no-install-recommends nettle-dev
RUN apt-get install -yq --no-install-recommends time
RUN apt-get install -yq --no-install-recommends libcairo2-dev
RUN apt-get install -yq --no-install-recommends gfortran
RUN apt-get install -yq --no-install-recommends libblas-dev
RUN apt-get install -yq --no-install-recommends libboost-all-dev
RUN apt-get install -yq --no-install-recommends libpstreams-dev
RUN apt-get install -yq --no-install-recommends cpanminus
RUN locale-gen en_US.UTF-8
RUN update-locale LANG=en_US.UTF-8

ENV OPT /opt/wtsi-cgp
ENV PATH $OPT/bin:$PATH
ENV PERL5LIB $OPT/lib/perl5
ENV LD_LIBRARY_PATH $OPT/lib
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV R_LIBS $OPT/R-libs
ENV R_LIBS_USER $R_LIBS
ENV R_PROFILE_USER $OPT/config/Rprofile
RUN mkdir -p $OPT
RUN mkdir $OPT/config/
RUN mkdir -p $R_LIBS_USER
WORKDIR /R
RUN cd /R && wget --progress=dot:giga -q -nv https://cdn.rstudio.com/r/ubuntu-2004/pkgs/r-4.2.1_1_amd64.deb &&\
	apt-get install -f -y ./r-4.2.1_1_amd64.deb &&\
	ln -s /opt/R/4.2.1/bin/R /usr/bin/R &&\
    ln -s /opt/R/4.2.1/bin/Rscript /usr/bin/Rscript &&\
    ln -s /opt/R/4.2.1/lib/R /usr/lib/R
WORKDIR /git
RUN cd /git && git clone --recursive https://github.com/cancerit/ascatNgs.git &&\
   cd ascatNgs &&\
   ./setup.sh $OPT

USER    ubuntu
FROM ubuntu:20.04

RUN apt-get -yq update && \
RUN apt-get install -qy --no-install-recommends \
lsb-release \
gnupg \
dirmngr \
wget \
build-essential \
locales \
gcc \
wget \
pkg-config \
cpanminus \
git \
unattended-upgrades && \
unattended-upgrade -d -v && \
apt-get remove -yq unattended-upgrades && \
apt-get autoremove -yq &&\
apt-get clean && \
rm -rf /var/lib/apt/lists/*
ENV DEBIAN_FRONTEND "noninteractive"
RUN locale-gen en_US.UTF-8
RUN update-locale LANG=en_US.UTF-8
RUN apt-get -yq update \
&& apt-get install -yq --no-install-recommends \
software-properties-common \
libreadline6-dev \
libx11-dev \
libxt-dev \
libpng-dev \
libjpeg-dev \
libcairo2-dev \
libcairo2-dev \
xvfb \
libbz2-dev \
libzstd-dev \
liblzma-dev \
libcurl4-openssl-dev \
texinfo \
texlive \
texlive-fonts-extra \
screen \
libpcre2-dev \
locales \
curl \
wget \
ca-certificates \
libperlio-gzip-perl \
bzip2 \
psmisc \
time \
zlib1g \
liblzma5 \
libncurses5 \
p11-kit \
libcurl3-gnutls \
libcurl4 \
moreutils \
google-perftools \
libcairo2 \
gfortran \
unattended-upgrades && \
unattended-upgrade -d -v && \
apt-get remove -yq unattended-upgrades && \
apt-get autoremove -yq && \
apt-get clean && \
rm -rf /var/lib/apt/lists/*

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN locale-gen en_US.UTF-8
RUN update-locale LANG=en_US.UTF-8

ENV OPT /opt/wtsi-cgp
ENV PATH $OPT/bin:$PATH
ENV PERL5LIB $OPT/lib/perl5
ENV LD_LIBRARY_PATH $OPT/lib
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV R_LIBS $OPT/R-libs
ENV R_LIBS_USER $R_LIBS
ENV R_PROFILE_USER $OPT/config/Rprofile

RUN mkdir -p $OPT
COPY --from=builder $OPT $OPT
COPY --from=builder /opt/R /opt/R
WORKDIR /R
RUN cd /R && wget --progress=dot:giga -q -nv https://cdn.rstudio.com/r/ubuntu-2004/pkgs/r-4.2.1_1_amd64.deb &&\
	apt-get install -f -y ./r-4.2.1_1_amd64.deb
RUN	ln -s /opt/R/4.2.1/bin/R /usr/bin/R
RUN ln -s /opt/R/4.2.1/bin/Rscript /usr/bin/Rscript
RUN ln -s /opt/R/4.2.1/lib/R /usr/lib/R
RUN rm -rf /R
RUN rm -rf /var/lib/apt/lists/*
RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu

USER    ubuntu
WORKDIR /home/ubuntu

CMD ["/bin/bash"]
