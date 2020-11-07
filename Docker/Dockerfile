# syntax = docker/dockerfile:1.0-experimental
FROM ubuntu:18.04 AS build
MAINTAINER Andrew Bailey, andbaile@ucsc.edu

# apt-get installs
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qq update && \
    apt-get -qq install -y --no-install-recommends autoconf wget git curl build-essential libbz2-dev zlib1g-dev liblzma-dev libeigen3-dev libreadline-gplv2-dev libncursesw5-dev libssl-dev libsqlite3-dev tk-dev libgdbm-dev libc6-dev libcurl4-openssl-dev ca-certificates python3.7-dev python3-pip python3.7-venv && \
    apt-get -qq clean && \
    apt-get -qq purge

# install cmake
WORKDIR /root
RUN mkdir /root/cmake && \
    wget -q https://cmake.org/files/v3.17/cmake-3.17.0-Linux-x86_64.sh --no-check-certificate && \
    sh /root/cmake-3.17.0-Linux-x86_64.sh --prefix=/root/cmake --skip-license && \
    ln -s /root/cmake/bin/cmake /usr/local/bin/cmake

# htslib
WORKDIR /root/
RUN wget -q https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 --no-check-certificate && \
    tar -vxjf htslib-1.9.tar.bz2 && \
    rm htslib-1.9.tar.bz2 && \
    cd /root/htslib-1.9 && \
    ./configure --prefix /usr/local --enable-plugins CPPFLAGS="-fPIC" CFLAGS="-fPIC" && \
    make && \
    make install

# boost install
WORKDIR /root/
RUN wget -q -O boost_1_69_0.tar.gz https://sourceforge.net/projects/boost/files/boost/1.69.0/boost_1_69_0.tar.gz/download --no-check-certificate && \
    tar -xzf boost_1_69_0.tar.gz >/dev/null && \
    rm boost_1_69_0.tar.gz && \
    cd /root/boost_1_69_0/ && \
    ./bootstrap.sh --with-libraries=system,date_time,filesystem,iostreams,coroutine,context,regex,thread,atomic >/dev/null && \
    ./b2 -d0 cxxflags="-fPIC" cflags="-fPIC" link=static -a && \
    ./b2 -d0 install

# hdf5 install
WORKDIR /root/
RUN wget -q https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.4/src/hdf5-1.10.4.tar.gz && \
    tar -xzf hdf5-1.10.4.tar.gz && \
    cd hdf5-1.10.4 && \
    ./configure --enable-threadsafe --disable-hl --prefix=/usr/local/ && \
    make && \
    make install

## setup python
WORKDIR /root/
RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && \
    python3.7 get-pip.py && \
    python3.7 -m pip -q install setuptools cython virtualenv


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

#RUN --mount=type=bind,source=.,target=/root/embed_fast5,rw cd /root/embed_fast5 && \
#    python3.7 -m pip install . && \
#    mkdir /opt/embed_fast5 && \
#    cp /root/embed_fast5/Docker/run_wrapper.sh /opt/embed_fast5/

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
#RUN --mount=type=bind,source=.,target=/root/embed_fast5,rw cd /root/embed_fast5 && \
#    pytest
COPY . /root/embed_fast5
RUN cd /root/embed_fast5 && \
    pytest
WORKDIR /data
