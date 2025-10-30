# Jasmine container installed via conda (bioconda)
# Lightweight base with micromamba for fast, reproducible installs
FROM mambaorg/micromamba:1.5.10

# Ensure commands run in bash and base env is active
SHELL ["/bin/bash", "-lc"]
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Install jasmine from bioconda (via conda-forge), then clean caches
RUN micromamba install -y -n base -c conda-forge -c bioconda jasminesv openjdk\
    && micromamba clean -a -y

# Optional: set a working directory
WORKDIR /data

# Expose conda binaries on PATH explicitly (usually already set)
ENV PATH=/opt/conda/bin:${PATH}

# Print version during build for quick verification (non-fatal)
RUN jasmine -h || true

# Default shell
CMD ["/bin/bash"]

