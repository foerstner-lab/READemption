From ubuntu:22.04

MAINTAINER Till Sauerwein <sauerwein@zbmed.de>

ARG DEBIAN_FRONTEND=noninteractive
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

ENV TZ=Europe/Berlin
RUN apt-get update && apt-get dist-upgrade --fix-missing -y
RUN apt-get -y install wget
RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda
RUN conda install -c till_sauerwein reademption -y

WORKDIR /root

