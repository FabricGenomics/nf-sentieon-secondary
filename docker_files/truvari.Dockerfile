FROM ubuntu:22.04

RUN apt-get -qq update \
  && DEBIAN_FRONTEND=noninteractive apt-get install -yq \
  curl \
  python3-dev \
  python3-pip \
  wget \
  bash \
  tabix \
  bcftools \
  && \
  rm -rf /var/lib/apt/lists/*

SHELL ["/bin/bash", "-c"]
WORKDIR /opt/truvari

RUN wget https://mafft.cbrc.jp/alignment/software/mafft_7.505-1_amd64.deb \
    && dpkg -i mafft_7.505-1_amd64.deb && rm mafft_7.505-1_amd64.deb

RUN python3 -m pip install --upgrade pip && \
    python3 -m pip install setproctitle pylint anybadge coverage && \
    python3 -m pip install --upgrade setuptools && \
    python3 -m pip install truvari==4.1.0

RUN mkdir -p /opt/truvari-source && \
    cd /opt/truvari-source && \
    wget -O adotto_TRregions_v1.2.1.bed.gz https://zenodo.org/records/8387564/files/adotto_TRregions_v1.2.bed.gz && \
    zcat adotto_TRregions_v1.2.1.bed.gz | cut -f1-3,18 | bgzip > anno.trf.bed.gz && \
    tabix anno.trf.bed.gz && \
    rm adotto_TRregions_v1.2.1.bed.gz && \
    wget -O trf409.linux64 https://github.com/Benson-Genomics-Lab/TRF/releases/download/v4.09.1/trf409.linux64 && \
    chmod +x trf409.linux64

WORKDIR /data
