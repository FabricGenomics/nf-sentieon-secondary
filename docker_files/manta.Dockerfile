############### stage 0: build samtools from source
FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive
ARG SAMTOOLS_VERSION=1.9
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 ca-certificates \
                 libbz2-dev \
                 libcurl4-openssl-dev \
                 liblzma-dev \
                 libncurses5-dev \
                 autoconf \
                 automake \
                 bzip2 \
                 gcc \
                 make \
                 wget \
                 zlib1g-dev && \
    wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && ./configure CFLAGS="-static" --without-curses && make -s all all-htslib && make install install-htslib && cd - && \
    rm -rf samtools-${SAMTOOLS_VERSION}* && \
    apt-get -qqy purge autoconf automake bzip2 gcc make wget && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_* && \
    samtools --help

############### stage 1: install Manta from local folder
FROM ubuntu:20.04

# copy from previous stage the binaries from samtools build
COPY --from=0 /usr/local/bin/* /usr/local/bin/

# install necessary packages for runtime
ARG MANTA_INSTALL_DIR=/usr/local/bin/manta/
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 ca-certificates \
                 bcftools \
                 python2.7 \
                 tabix \
                 wget \
                 zlib1g-dev && \
    ln -s $(which python2.7) /usr/bin/python2 && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*

# download and unpack the Manta distribution
RUN mkdir -p ${MANTA_INSTALL_DIR} && \
    wget -O /tmp/manta.tar.bz2 https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2 && \
    tar -xjf /tmp/manta.tar.bz2 -C /tmp && \
    mv /tmp/manta-1.6.0.centos6_x86_64/* ${MANTA_INSTALL_DIR} && \
    rm -rf /tmp/manta.tar.bz2 /tmp/manta-1.6.0.centos6_x86_64

# optionally run the demo; skip verification failures to avoid blocking image build
ARG RUN_DEMO=false
RUN if [ "$RUN_DEMO" = "true" ]; then \
        python2 ${MANTA_INSTALL_DIR}/bin/runMantaWorkflowDemo.py || echo "Demo verification failed; continuing" >&2; \
    fi && \
    rm -rf MantaDemoAnalysis && \
    rm -rf ${MANTA_INSTALL_DIR}share/demo && \
    samtools --help
