FROM ubuntu:16.04
# From https://github.com/biotoolsdx/docker-gffcompare-gffread/blob/master/Dockerfile

RUN apt-get update && apt-get install -y build-essential cmake git && \
    mkdir -p /opt/

WORKDIR /opt/
RUN git clone https://github.com/gpertea/gclib \
    && git clone https://github.com/gpertea/gffcompare \
    && git clone https://github.com/gpertea/gffread \
    && cd gffcompare \
    && make release \
    && cd /opt/gffread \
    && make release

ENV PATH "/opt/gffcompare:/opt/gffread:$PATH"

CMD bash