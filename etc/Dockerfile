FROM ubuntu:24.04 as base
LABEL authors="nicholashathaway"

ARG DEBIAN_FRONTEND="noninteractive"
ARG SHELL="/bin/bash"
ARG LANG="en_US.UTF-8"
ARG LANGUAGE="en_US.UTF-8"
ARG LC_ALL="en_US.UTF-8"
ARG CPU_COUNT=15
ARG TIME_ZONE=Etc/UTC
ARG ELUCIDATOR_VERSION=develop

RUN ulimit -n 10000

# install wget for getting pubkey
RUN apt-get update
RUN apt-get -yq dist-upgrade wget

# add packages source to install r 4.4 which is needed for BiocManager 3.19 which is needed pwalign
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu noble-cran40/" > /etc/apt/sources.list.d/cran.list
RUN echo "deb-src https://cloud.r-project.org/bin/linux/ubuntu noble-cran40/" >> /etc/apt/sources.list.d/cran.list
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc

# install basics
RUN apt-get update
RUN apt-get -yq dist-upgrade
RUN apt-get install -yq --no-install-recommends autotools-dev autoconf libtool automake build-essential curl wget file tree\
    git locales libssl-dev libcurl4-gnutls-dev ca-certificates xz-utils zlib1g-dev libbz2-dev liblzma5 liblzma-doc\
    liblzma-dev openssh-client gcc-10 g++-10 ccache git-lfs libxml2-dev libfontconfig1-dev libcairo2-dev\
    python3-dev python3-pip python3-nbclient python3-nbconvert python3-ipykernel python3-jupyter-client python3-jupyter-core

# install common bio tools
RUN apt-get install -yq --no-install-recommends lastz ncbi-blast+ ncbi-tools-bin hmmer minimap2 bwa bowtie2 cmake  \
    samtools mummer muscle pigz r-base r-base-dev

# clean up 
RUN apt-get -yq autoclean
RUN apt-get -yq autoremove

# install snakemake
RUN pip3 install --break-system-packages snakemake

# set environment locale
RUN echo "$LANG UTF-8" >> /etc/locale.gen
RUN echo "LANG=$LANG" > /etc/locale.conf
RUN echo "LC_ALL=$LC_ALL" >> /etc/environment
RUN echo "LANGUAGE=$LANGUAGE" >> /etc/environment
RUN locale-gen $LANG
RUN update-locale LANG=$LANG

RUN DEBIAN_FRONTEND=noninteractive TZ=$TIME_ZONE apt-get -y install tzdata

# add github key so can git clone from github
RUN mkdir ~/.ssh/
RUN ssh-keyscan github.com >> ~/.ssh/known_hosts

# install bwa-mem2
WORKDIR /opt
RUN curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 | tar jxf -

# install quarto
WORKDIR /opt
RUN wget -q https://github.com/quarto-dev/quarto-cli/releases/download/v1.6.40/quarto-1.6.40-linux-amd64.tar.gz
RUN tar -xzf quarto-1.6.40-linux-amd64.tar.gz
RUN rm quarto-1.6.40-linux-amd64.tar.gz
RUN /opt/quarto-1.6.40/bin/quarto check

# pmotools
WORKDIR /opt
RUN git clone  https://github.com/PlasmoGenEpi/pmotools-python.git
WORKDIR /opt/pmotools-python
RUN git checkout develop
RUN pip install --break-system-packages .

# install R packages
WORKDIR /opt
RUN mkdir -p /usr/local/lib/R/etc/ /usr/lib/R/etc/
RUN echo "options(repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl',  \
    Ncpus = ${CPU_COUNT})" | tee /usr/local/lib/R/etc/Rprofile.site | tee /usr/lib/R/etc/Rprofile.site
RUN R -e 'install.packages(c("remotes"))'
## attempt to load libraries to make sure they installed
RUN R -e 'library("remotes")'


# R packages
RUN Rscript -e "install.packages(c( \
     'tibble', 'dplyr', 'stringr', 'readr', 'ggplot2', 'tidyr', \
     'data.table', 'openxlsx', 'rlang', 'magrittr', 'pegas',  \
     'ape', 'rmarkdown', 'knitr', 'kableExtra', 'DT','here', \
     'ggrepel', 'ggthemes', 'ggdist', 'ggnewscale', 'plotly', \
     'dendsort', 'tidymodels', 'rjson', 'cowplot', 'dbscan', 'naturalsort', 'Cairo', \
     'fastmatch', 'colorRamps', 'paletteer','r2d3', 'shiny', 'skimr', 'moderndive', \
     'formatR', 'markdown', 'bigmds'), Ncpus = ${CPU_COUNT}, repos = 'http://cran.us.r-project.org')"

## attempt to load libraries to make sure they installed
RUN R -e 'library("tibble")'
RUN R -e 'library("dplyr")'
RUN R -e 'library("stringr")'
RUN R -e 'library("readr")'
RUN R -e 'library("ggplot2")'
RUN R -e 'library("tidyr")'
RUN R -e 'library("data.table")'
RUN R -e 'library("openxlsx")'
RUN R -e 'library("rlang")'
RUN R -e 'library("magrittr")'
RUN R -e 'library("pegas")'
RUN R -e 'library("ape")'
RUN R -e 'library("rmarkdown")'
RUN R -e 'library("knitr")'
RUN R -e 'library("kableExtra")'
RUN R -e 'library("DT")'
RUN R -e 'library("here")'
RUN R -e 'library("ggrepel")'
RUN R -e 'library("ggthemes")'
RUN R -e 'library("ggdist")'
RUN R -e 'library("ggnewscale")'
RUN R -e 'library("plotly")'
RUN R -e 'library("dendsort")'
RUN R -e 'library("tidymodels")'
RUN R -e 'library("rjson")'
RUN R -e 'library("cowplot")'
RUN R -e 'library("dbscan")'
RUN R -e 'library("naturalsort")'
RUN R -e 'library("Cairo")'
RUN R -e 'library("fastmatch")'
RUN R -e 'library("colorRamps")'
RUN R -e 'library("paletteer")'
RUN R -e 'library("r2d3")'
RUN R -e 'library("shiny")'
RUN R -e 'library("skimr")'
RUN R -e 'library("moderndive")'
RUN R -e 'library("formatR")'
RUN R -e 'library("markdown")'
RUN R -e 'library("bigmds")'

# Bioconductor packages (have to install 3.19 because that's needed for pwalign)
RUN Rscript -e 'if (!require("BiocManager", quietly = TRUE)) { install.packages("BiocManager"); };  \
    BiocManager::install(version = "3.19", ask = FALSE);'
RUN Rscript -e 'BiocManager::install("Biostrings", ask = FALSE)'
RUN Rscript -e 'BiocManager::install("pwalign", ask = FALSE)'
RUN Rscript -e 'BiocManager::install("ggtree", ask = FALSE)'
RUN Rscript -e 'BiocManager::install("ggtreeExtra", ask = FALSE)'
RUN Rscript -e 'BiocManager::install("ComplexHeatmap", ask = FALSE)'

## attempt to load libraries to make sure they installed
RUN R -e 'library("Biostrings")'
RUN R -e 'library("pwalign")'
RUN R -e 'library("ggtree")'
RUN R -e 'library("ggtreeExtra")'
RUN R -e 'library("ComplexHeatmap")'

# install github r packages
RUN Rscript -e "remotes::install_github('nickjhathaway/HaplotypeRainbows@develop')"
RUN Rscript -e "remotes::install_github('hoesler/rwantshue')"
RUN R -e 'library("HaplotypeRainbows")'
RUN R -e 'library("rwantshue")'

# install elucidator
WORKDIR /opt
RUN git clone https://github.com/nickjhathaway/elucidator.git
WORKDIR /opt/elucidator
RUN git checkout develop
RUN ./setup.py  --libs adapterremoval:v2.3.2 --numCores ${CPU_COUNT} --overWrite --symlinkBin
RUN ./setup.py  --libs bamtools:develop --symlinkBin
RUN ./setup.py  --libs seqServer:develop --numCores ${CPU_COUNT} --overWrite --symlinkBin
RUN ./setup.py  --libs SeekDeep:develop --numCores ${CPU_COUNT} --overWrite --symlinkBin
RUN ./setup.py  --libs PathWeaver:develop --numCores ${CPU_COUNT} --overWrite --symlinkBin
RUN ./setup.py  --libs elucidator:${ELUCIDATOR_VERSION} --numCores ${CPU_COUNT} --overWrite --symlinkBin
RUN rm -fr external/build

# expose ports
EXPOSE 8000-9999

# update path
ENV PATH "/opt/pmotools-python/scripts:/opt/quarto-1.6.40/bin/:/opt/bwa-mem2-2.2.1_x64-linux:/opt/elucidator/bin:$PATH"
ENV XDG_CACHE_HOME "/tmp/quarto_cache"
# reset working directory
WORKDIR /
