From ubuntu:20.04

MAINTAINER Till Sauerwein <sauerwein@zbmed.de>

ARG DEBIAN_FRONTEND=noninteractive
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

ENV TZ=Europe/Berlin
RUN apt-get update && apt-get dist-upgrade -y
RUN apt-get install -y python3 python3-setuptools python3-pip python3-matplotlib cython3 zlib1g-dev make libncurses5-dev r-base libxml2-dev
RUN apt-get install wget tree
RUN pip install READemption
RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh
RUN conda --version
RUN conda update --yes conda
RUN conda install -c bioconda segemehl=0.3.4
RUN echo "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager') ; BiocManager::install('DESeq2'); install.packages('gplots')" | R --no-save

WORKDIR /root

