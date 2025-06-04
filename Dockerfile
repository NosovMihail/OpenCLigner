FROM ubuntu:20.04

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    build-essential \
    ocl-icd-opencl-dev \
    clinfo \
    cmake \
    git \
    nano

WORKDIR /app
COPY . /app

RUN g++ -std=c++17 -O2 -o opencligner main.cpp -lOpenCL