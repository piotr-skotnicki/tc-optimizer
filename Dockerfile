# TC Optimizing Compiler

FROM ubuntu:22.04

LABEL maintainer="pskotnicki@zut.edu.pl"

RUN apt-get update && apt-get install -y automake autoconf libtool pkg-config libgmp3-dev \
                        libclang-dev llvm libntl-dev g++ make git clang zlib1g-dev libglpk-dev

WORKDIR /app

RUN git clone https://github.com/piotr-skotnicki/tc-optimizer.git tc

WORKDIR /app/tc

RUN git submodule update --init --recursive \
    && ./autogen.sh \
    && ./configure \
    && make

WORKDIR /app/tc/src

CMD ["./tc", "--help"]
