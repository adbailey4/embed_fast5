# syntax = docker/dockerfile:1.0-experimental
FROM baileybase/embed_dependencies:v0.0.2 AS build
MAINTAINER Andrew Bailey, andbaile@ucsc.edu

ENV VIRTUAL_ENV=/opt/venv
RUN python3.7 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"
# build static so we need to tranfser less shared dependencies
RUN python3.7 -m pip -q install wheel

ENV BUILD_SHARED_LIBS="OFF"
COPY . /root/embed_fast5
RUN cd /root/embed_fast5 && \
    python3.7 -m pip install . && \
    mkdir /opt/embed_fast5 && \
    cp /root/embed_fast5/Docker/run_wrapper.sh /opt/embed_fast5/

FROM ubuntu:18.04 as runtime
RUN apt-get -qq update && \
    apt-get -qq install -y --no-install-recommends python3.7 python3.7-venv libgomp1 && \
    apt-get -qq clean && \
    apt-get -qq purge
WORKDIR /home/

COPY --from=build /opt/embed_fast5 /opt/embed_fast5/
COPY --from=build /opt/venv /opt/venv

ENV VIRTUAL_ENV=/opt/venv
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

RUN mkdir /data
WORKDIR /data

ENTRYPOINT ["sh", "/opt/embed_fast5/run_wrapper.sh"]

FROM runtime as test
RUN apt-get -qq install -y --no-install-recommends python3-distutils
RUN python3.7 -m pip install pytest
COPY . /root/embed_fast5
RUN cd /root/embed_fast5 && \
    pytest
WORKDIR /data
