# TC Optimizing Compiler

FROM ubuntu:20.04

LABEL maintainer="pskotnicki@zut.edu.pl"

RUN apt-get update
RUN apt-get install -y automake autoconf libtool pkg-config libgmp3-dev libclang-dev \
                       llvm libntl-dev g++ make git clang zlib1g-dev libglpk-dev
RUN cd /tmp \
    && git clone https://github.com/piotr-skotnicki/tc-optimizer.git tc \
    && cd tc \
    && git submodule update --init --recursive \
    && ./autogen.sh \
    && ./configure \
    && make
