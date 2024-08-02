FROM ubuntu:22.04 as base
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

# install basics
RUN apt-get update
RUN apt-get -yq dist-upgrade
RUN apt-get install -yq --no-install-recommends autotools-dev autoconf libtool automake build-essential curl wget file git locales libssl-dev libcurl4-gnutls-dev ca-certificates xz-utils zlib1g-dev libbz2-dev liblzma5 liblzma-doc liblzma-dev openssh-client gcc-10 g++-10 python3-dev python3-pip ccache

# install common bio tools
RUN apt-get install -yq --no-install-recommends lastz ncbi-blast+ ncbi-tools-bin hmmer minimap2 bwa bowtie2 cmake samtools mummer muscle

# install snakemake
RUN pip3 install snakemake

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

WORKDIR /opt
RUN git clone https://github.com/nickjhathaway/elucidator.git
WORKDIR /opt/elucidator
RUN git checkout develop
RUN ./setup.py  --libs adapterremoval:v2.3.2 --numCores ${CPU_COUNT} --overWrite --symlinkBin
RUN ./setup.py  --libs elucidator:${ELUCIDATOR_VERSION} --numCores ${CPU_COUNT} --overWrite --symlinkBin

# expose ports
EXPOSE 8000-9999


# update path
ENV PATH "/opt/elucidator/bin:$PATH"

# reset working directory
WORKDIR /


