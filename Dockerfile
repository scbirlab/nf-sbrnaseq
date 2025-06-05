FROM continuumio/miniconda3

COPY environment.yml /tmp/environment.yml

RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    make \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

RUN conda env create -n env -f /tmp/environment.yml && \
    conda clean --all --yes

SHELL ["conda", "run", "-n", "env", "/bin/bash", "-c"]
ENV PATH=/opt/conda/envs/env/bin:$PATH
