#
# docker build -t oldnoselecs .
#
# docker run --name oldnoselecs -it --rm -v ~/Work/Circuit/ubuntu_10.04/src:/src oldnoselecs /bin/bash
#
# su pissard
# 


FROM ubuntu:10.04
RUN sed -i 's/archive/old-releases/g' /etc/apt/sources.list
RUN apt-get -y update
RUN apt-get -y install emacs vim git-core
RUN apt-get -y install cmake  gcc gfortran g++ liblapack-dev libboost-dev libboost-graph-dev   libxml2-dev libgmp3-dev
RUN apt-get -y install libx11-dev libxt-dev
VOLUME /src


RUN addgroup --gid 50101 sed\
    && useradd pissard --uid 50107 --gid 50101 -d /home/pissard \
    && mkdir /home/pissard \
    && chown pissard:sed /home/pissard

RUN addgroup --gid 50096 bipop\
    && useradd acary --uid 51029  --gid 50096 -d /home/pissard \
    && mkdir /home/acary \
    && chown acary:bipop /home/acary
   
WORKDIR /src
