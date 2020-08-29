FROM gcc:6 AS build

MAINTAINER Andrew Bailey, andbaile@ucsc.edu

# apt-get installs
RUN apt-get update -qq && apt-get install -y --no-install-recommends build-essential libbz2-dev zlib1g-dev liblzma-dev libeigen3-dev
# cmake
WORKDIR /root/
RUN wget https://github.com/Kitware/CMake/releases/download/v3.17.4/cmake-3.17.4.tar.gz && tar -xzvf cmake-3.17.4.tar.gz && rm cmake-3.17.4.tar.gz
WORKDIR /root/cmake-3.17.4
RUN ./bootstrap && make && make install

# htslib
WORKDIR /root/
RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && tar -vxjf htslib-1.9.tar.bz2 && rm htslib-1.9.tar.bz2
WORKDIR /root/htslib-1.9
RUN ./configure --prefix /usr/local --enable-plugins CPPFLAGS="-fPIC" CFLAGS="-fPIC" && make && make install

# boost install
ENV CFLAGS "$CFLAGS -fPIC"
ENV CXXFLAGS "$CXXFLAGS -fPIC"
WORKDIR /root/
RUN wget -O boost_1_69_0.tar.gz https://sourceforge.net/projects/boost/files/boost/1.69.0/boost_1_69_0.tar.gz/download && tar -xzf boost_1_69_0.tar.gz >/dev/null && rm boost_1_69_0.tar.gz
WORKDIR /root/boost_1_69_0/
#RUN ./bootstrap.sh --show-libraries
RUN ./bootstrap.sh --with-libraries=system,date_time,filesystem,iostreams,coroutine,context >/dev/null && ./b2 && ./b2 install

# hdf5
WORKDIR /root/
RUN wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.4/src/hdf5-1.10.4.tar.gz && tar -xzf hdf5-1.10.4.tar.gz && rm hdf5-1.10.4.tar.gz
RUN cd hdf5-1.10.4 && ./configure --prefix /usr/local --enable-threadsafe --disable-hl && make && make install
# python 3 install
WORKDIR /root/
RUN wget https://www.python.org/ftp/python/3.6.8/Python-3.6.8.tgz && tar -xzf Python-3.6.8.tgz && rm Python-3.6.8.tgz
WORKDIR /root/Python-3.6.8/
RUN ./configure && make && make install
## cython
RUN pip3 install cython

#/root
## embed
WORKDIR /root/
#RUN ls
#RUN git clone --recursive https://github.com/adbailey4/embed_fast5.git --branch filter_by_position
COPY . /root/embed_fast5
RUN mkdir build2
WORKDIR /root/embed_fast5/build2

#WORKDIR /embed_fast5

#RUN mkdir build
#WORKDIR /home/embed_fast5/build
#RUN cmake ..

#RUN python3 setup.py install
#
#ENV LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib"
#
#RUN embed_main
#RUN pip3 install ecdsa==0.15
#RUN python3 setup.py test
#RUN embed_main
#
#FROM alpine:latest as runtime
#RUN apk update && apk add --no-cache \
#    libstdc++
#WORKDIR /home/
#COPY --from=build /usr/local/bin/embed_main /bin/embed_main
##COPY --from=build /usr/ /usr/
#RUN /bin/embed_main
#
## set signalAlign bin as workDir
##WORKDIR /home/signalAlign/bin/
##
#COPY run_wrapper.sh /opt/embed_fast5/
#
#RUN mkdir /data
#WORKDIR /data
#
#ENTRYPOINT ["sh", "/opt/embed_fast5/run_wrapper.sh"]
