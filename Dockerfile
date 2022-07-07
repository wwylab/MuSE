FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y git g++ cmake autoconf libtool liblzma-dev zlib1g-dev libbz2-dev libcurl3-dev libssl-dev

ADD . /MuSE
# Remove CR from DOS if running in Windows
RUN cd MuSE && tr -d '\r' < install_muse.sh > install_muse_unix.sh
RUN cd MuSE && bash ./install_muse_unix.sh

RUN mkdir /MuSE/bin
RUN cp /MuSE/MuSE /MuSE/bin
ENV PATH=$PATH:/MuSE/bin
