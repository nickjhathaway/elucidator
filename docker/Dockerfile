FROM ubuntu:24.04 AS base
LABEL authors="nicholashathaway"

ARG DEBIAN_FRONTEND="noninteractive"
ARG CPU_COUNT=15
ARG TIME_ZONE=Etc/UTC
ARG ELUCIDATOR_VERSION=develop
ARG LANG="en_US.UTF-8"

RUN ulimit -n 10000


# All apt installs in one layer — install everything upfront, clean once
RUN apt-get update && apt-get -yq dist-upgrade && \
    apt-get install -yq --no-install-recommends \
      autotools-dev autoconf libtool automake build-essential curl wget file tree \
      git locales libssl-dev libcurl4-gnutls-dev ca-certificates \
      xz-utils zlib1g-dev libbz2-dev liblzma5 liblzma-dev openssh-client \
      gcc-10 g++-10 ccache libxml2-dev libfontconfig1-dev libcairo2-dev \
      python3-dev python3-pip python3-yaml \
      libharfbuzz-dev libfribidi-dev libicu-dev libhdf5-dev \
      pigz lastz cmake tzdata libuv1-dev\
      # R and CRAN repo
      r-base r-base-dev \
 && echo "deb https://cloud.r-project.org/bin/linux/ubuntu noble-cran40/" > /etc/apt/sources.list.d/cran.list \
 && wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
      | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
 && apt-get update && apt-get install -yq r-base r-base-dev \
 && pip3 install --break-system-packages --no-cache-dir --ignore-installed setuptools wheel \
 && rm -rf /usr/lib/python3/dist-packages/setuptools* /usr/lib/python3/dist-packages/wheel* \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

# Locale + timezone in one layer
RUN echo "$LANG UTF-8" >> /etc/locale.gen && \
    locale-gen $LANG && update-locale LANG=$LANG && \
    ln -fs /usr/share/zoneinfo/$TIME_ZONE /etc/localtime && \
    dpkg-reconfigure -f noninteractive tzdata

# Miniforge + all conda tools in one layer
RUN wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \
      -O /tmp/miniforge.sh \
 && bash /tmp/miniforge.sh -b -p /opt/conda \
 && rm /tmp/miniforge.sh \
 && /opt/conda/bin/mamba install -y -c bioconda -c conda-forge \
      muscle minimap2 bowtie2 bwa samtools bcftools \
      bedtools seqkit hmmer blast salmon ncbi-datasets-cli \
      gffread subread mosdepth fastp seqtk mummer \
 && /opt/conda/bin/mamba clean -afy \
 && find /opt/conda -name "*.pyc" -delete \
 && find /opt/conda -name "*.a" -delete \
 # drop build-time artifacts we never use at runtime (only prebuilt tools run)
 && rm -rf /opt/conda/x86_64-conda-linux-gnu \
 && rm -rf /opt/conda/include \
 && rm -rf /opt/conda/pkgs/* /opt/conda/share/doc /opt/conda/share/man /opt/conda/share/info \
 # drop conda's python + conda/mamba machinery — unused at runtime: no compiled tool
 # links libpython (verified), and the env-python3 helper scripts (bowtie2/blast) fall
 # back to the system python 3.12.
 && rm -rf /opt/conda/lib/python3.* /opt/conda/lib/libpython3.* \
           /opt/conda/bin/python* /opt/conda/bin/pip* \
           /opt/conda/bin/conda* /opt/conda/bin/mamba* /opt/conda/condabin

ENV PATH="/opt/conda/bin:$PATH"

# All R packages in one layer — only what the etc/*.qmd reports require
COPY install_packages.R /opt/install_packages.R
RUN mkdir -p /usr/local/lib/R/etc/ /usr/lib/R/etc/ && \
    echo "options(repos = c(CRAN = 'https://cran.rstudio.com/'), \
      download.file.method = 'libcurl', Ncpus = ${CPU_COUNT})" \
      | tee /usr/local/lib/R/etc/Rprofile.site | tee /usr/lib/R/etc/Rprofile.site && \
    Rscript -e " \
      source('/opt/install_packages.R'); \
      install_cran(c( \
        'remotes', 'readr', 'ggplot2', 'dplyr', 'tidyr', 'stringr', \
        'DT', 'plotly', 'rmarkdown', 'knitr', 'Cairo' \
      ), ncpus = ${CPU_COUNT})" && \
    Rscript -e " \
      if (!require('BiocManager', quietly=TRUE)) install.packages('BiocManager'); \
      BiocManager::install(version='3.23', ask=FALSE);" && \
    Rscript -e "\
      source('/opt/install_packages.R'); \
      install_bioconductor(c('ComplexHeatmap'), version = '3.23')" && \
    Rscript -e " \
      source('/opt/install_packages.R'); \
      setRepositories(ind = 1:3);\
      install_github(c('nickjhathaway/HaplotypeRainbows@develop', 'hoesler/rwantshue'))" && \
    # Clean R build artifacts — often 200-400MB
    find /usr/local/lib/R -name "*.o" -delete && \
    find /usr/local/lib/R -name "*.a" -delete && \
    rm -rf /tmp/downloaded_packages /tmp/Rtmp* \
    /root/.local/share/R



# Quarto
WORKDIR /opt
RUN wget -q https://github.com/quarto-dev/quarto-cli/releases/download/v1.9.38/quarto-1.9.38-linux-amd64.tar.gz && \
    tar -xzf quarto-1.9.38-linux-amd64.tar.gz && \
    rm quarto-1.9.38-linux-amd64.tar.gz && \
    # quarto's bundled esbuild is built with an old Go toolchain (Scout stdlib CVEs). esbuild is
    # still required (dashboard/embed-resources rendering spawns it to bundle JS), so swap in the
    # latest upstream esbuild binary, which is built with a current, CVE-free Go stdlib.
    wget -q https://registry.npmjs.org/@esbuild/linux-x64/-/linux-x64-0.28.1.tgz -O /tmp/esbuild.tgz && \
    tar -xzOf /tmp/esbuild.tgz package/bin/esbuild > /opt/quarto-1.9.38/bin/tools/x86_64/esbuild && \
    chmod +x /opt/quarto-1.9.38/bin/tools/x86_64/esbuild && \
    rm /tmp/esbuild.tgz

# elucidator — shallow clone, clean build artifacts
RUN mkdir ~/.ssh/ && ssh-keyscan github.com >> ~/.ssh/known_hosts
RUN git clone --depth 1 -b ${ELUCIDATOR_VERSION} https://github.com/nickjhathaway/elucidator.git /opt/elucidator
WORKDIR /opt/elucidator
RUN ./setup.py --libs adapterremoval:v2.3.2 --numCores ${CPU_COUNT} --overWrite --symlinkBin && \
    # elucidator's setup.py resolves and builds the correct pinned versions of its own
    # deps (bamtools, njhcpp, njhseq, seqServer, SeekDeep, PathWeaver, ...), so installing
    # them explicitly here would clone/build the wrong (develop) versions for a tagged build
    ./setup.py --libs elucidator:${ELUCIDATOR_VERSION} --numCores ${CPU_COUNT} --overWrite && \
    # --symlinkBin only links the single tool it's building; instead link every built tool's
    # binaries (elucidator, SeekDeep, PathWeaver, deps) into one dir on PATH. All tools install
    # their binaries under external/local/<tool>/<version>/<build>/bin/
    mkdir -p /opt/elucidator/bin && \
    find -L /opt/elucidator/external/local/*/*/*/bin/ -maxdepth 1 -type f -executable \
      -exec ln -sft /opt/elucidator/bin/ {} + && \
    rm -fr external/build external/src && \
    find . -name "*.o" -delete && \
    find /opt/elucidator/external/local -name ".git" -type d  | xargs rm -rf && \
    find /opt/elucidator/external/local -name "tests" -type d  | xargs rm -rf && \
    find /opt/elucidator/external/local -name "docs" -type d  | xargs rm -rf && \
    find /opt/elucidator/external/local -name "doc" -type d  | xargs rm -rf


EXPOSE 8000-9999
ENV PATH="/opt/quarto-1.9.38/bin:/opt/elucidator/bin:$PATH"
ENV XDG_CACHE_HOME="/tmp/quarto_cache"
WORKDIR /

# Default to a non-root user (Docker Scout policy). This is only a default for a bare `docker run`:
# Nextflow's `-u $(id -u):$(id -g)` overrides it at runtime (so file ownership is unchanged), and
# Apptainer/Singularity ignore USER entirely and run as the calling user. Tools live in world-readable
# /opt and the Quarto cache is under world-writable /tmp, so any uid can run the image.
RUN useradd --create-home --uid 1000 --shell /bin/bash elucidator
USER elucidator
