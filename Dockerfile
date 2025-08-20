FROM mambaorg/micromamba:1.5.6

# Add a dummy user entry to avoid `passwd` errors in Docker
USER root
RUN echo "user:x:1001:1001::/home/user:/bin/bash" >> /etc/passwd && \
    mkdir -p /home/user && chown -R 1001:1001 /home/user

COPY environment.yml /tmp/environment.yml
RUN micromamba create -n env -f /tmp/environment.yml && \
    micromamba clean --all --yes
RUN mkdir -p /tmp/{xdg,fontconfig,mpl,numba}

ENV PATH=/opt/conda/envs/env/bin:$PATH
ENV MPLBACKEND=Agg
ENV XDG_CACHE_HOME=/tmp/xdg
ENV MPLCONFIGDIR=/tmp/mpl
ENV NUMBA_CACHE_DIR=/tmp/numba
