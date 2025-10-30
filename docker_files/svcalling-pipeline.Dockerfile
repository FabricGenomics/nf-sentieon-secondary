# syntax=docker/dockerfile:1.4

# Build stage: compile PanGenie from local source
FROM --platform=linux/amd64 ubuntu:24.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive \
    LC_ALL=C

ARG PANGENIE_REPO=https://github.com/eblerjana/pangenie.git
ARG PANGENIE_REF=master

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        git \
        build-essential \
        zlib1g-dev \
        libcereal-dev \
        libjellyfish-2.0-dev \
        pkg-config \
        cmake \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
RUN git clone --branch "${PANGENIE_REF}" --depth 1 "${PANGENIE_REPO}" /src/pangenie

# Build PanGenie
WORKDIR /src/pangenie
RUN mkdir -p build \
    && cd build \
    && cmake .. \
    && make -j"$(nproc)"

# Collect binaries and metadata for the runtime image
RUN mkdir -p /out/bin /out/metadata \
    && cp -v build/src/PanGenie /out/bin/ \
    && cp -v build/src/PanGenie-index /out/bin/ \
    && cp -v build/src/PanGenie-vcf /out/bin/ \
    && cp -v build/src/PanGenie-sampling /out/bin/ \
    && cp -v build/src/Analyze-UK /out/bin/ \
    && mkdir -p /out/lib \
    && cp -v build/src/libPanGenieLib.so /out/lib/ \
    && bash -lc 'dpkg -l | grep jellyfish | tr -s " " | cut -d " " -f 2,3 > /out/metadata/jellyfish.lib.version || true' \
    && bash -lc 'git -C /src/pangenie rev-parse --short HEAD > /out/metadata/pangenie.git.version || echo unknown > /out/metadata/pangenie.git.version'


# Runtime stage: minimal deps to run the binaries
FROM --platform=linux/amd64 ubuntu:24.04

ENV LC_ALL=C \
    PANGENIE_HOME=/repos \
    LD_LIBRARY_PATH=/usr/local/lib:${LD_LIBRARY_PATH}

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        zlib1g \
        python3 \
        python3-pip \
        python3-pysam \
        python3-requests \
        python3-boto3 \
        perl \
        libjson-perl \
        libjellyfish-2.0-2 \
        tabix \
        samtools \
        bcftools \
    && rm -rf /var/lib/apt/lists/*

# Add compiled binaries and metadata
COPY --from=builder /out/bin/ /usr/local/bin/
COPY --from=builder /out/lib/ /usr/local/lib/
COPY --from=builder /out/metadata/ /metadata/

# Include the PanGenie repository cloned during build for reference
RUN mkdir -p ${PANGENIE_HOME}
COPY --from=builder /src/pangenie ${PANGENIE_HOME}/pangenie
RUN rm -rf ${PANGENIE_HOME}/pangenie/build

