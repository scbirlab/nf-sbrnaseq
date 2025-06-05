FROM mambaorg/micromamba:1.5.6

COPY environment.yml /tmp/environment.yml

# Install build tools required for pip installs
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    make \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Create conda environment
RUN micromamba create -n env -f /tmp/environment.yml && \
    micromamba clean --all --yes

ENV PATH=/opt/conda/envs/env/bin:$PATH
