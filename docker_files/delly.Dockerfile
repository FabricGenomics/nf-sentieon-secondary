# use the ubuntu base image
FROM ubuntu:22.04

LABEL maintainer="Tobias Rausch <rausch@embl.de>"

# install required packages using apt
RUN apt-get update && apt-get install -y \
    autoconf \
    automake \
    build-essential \
    cmake \
    g++ \
    gfortran \
    git \
    libcurl4-gnutls-dev \
    libtool \
    hdf5-tools \
    libboost-date-time-dev \
    libboost-program-options-dev \
    libboost-system-dev \
    libboost-filesystem-dev \
    libboost-iostreams-dev \
    libbz2-dev \
    libdeflate-dev \
    libhdf5-dev \
    libncurses-dev \
    liblzma-dev \
    pkg-config \
    zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# set environment
ENV BOOST_ROOT=/usr

# copy local delly repo and build
# install delly from git (recursive to include htslib)
RUN cd /opt \
    && git clone --recursive https://github.com/dellytools/delly.git \
    && cd /opt/delly/ \
    && make STATIC=1 all \
    && make install


# Multi-stage build
FROM ubuntu:22.04
RUN mkdir -p /opt/delly/bin /maps
WORKDIR /opt/delly/bin
COPY --from=0 /opt/delly/bin/delly .

# Workdir
WORKDIR /home

# Add Delly to PATH
ENV PATH="/opt/delly/bin:${PATH}"

# Get bcftools and tabix
RUN apt-get update \
  && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
     bcftools \
     tabix \
     wget \
     ca-certificates \
  && rm -rf /var/lib/apt/lists/*

# Copy reference files into /maps inside the image
RUN cd /maps && \
    wget https://gear-genomics.embl.de/data/delly/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz && \
    wget https://gear-genomics.embl.de/data/delly/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz.fai && \
    wget https://gear-genomics.embl.de/data/delly/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz.gzi
