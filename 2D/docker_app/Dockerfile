FROM ubuntu

RUN apt update
RUN apt install -y gcc-13-offload-nvptx
RUN apt install -y g++
RUN apt install -y make

WORKDIR /rp
RUN /bin/bash