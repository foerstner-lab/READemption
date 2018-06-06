From ubuntu:16.04
MAINTAINER Konrad Förstner <konrad.foerstner@uni-wuerzburg.de>

RUN apt-get update && apt-get dist-upgrade -y
RUN apt-get install -y python3 python3-setuptools python3-pip python3-matplotlib cython3 zlib1g-dev make libncurses5-dev r-base libxml2-dev

RUN curl http://www.bioinf.uni-leipzig.de/Software/segemehl/segemehl_0_2_0.tar.gz > segemehl_0_2_0.tar.gz && \
    tar xzf segemehl_0_2_0.tar.gz && \
    cd segemehl_*/segemehl/ && make && cd ../../ && \
    cp segemehl_*/segemehl/segemehl.x /usr/bin/ && \
    cp segemehl_*/segemehl/lack.x /usr/bin/  && \
    rm -rf segemehl*

RUN echo 'source("http://bioconductor.org/biocLite.R"); biocLite("DESeq2")' | R --no-save
RUN pip3 install READemption

WORKDIR /root

