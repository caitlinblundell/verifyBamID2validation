# Define base image Ubuntu 20.04 as Build Stage
FROM ubuntu:20.04 AS build

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive
ENV LD_LIBRARY_PATH=/usr/local/lib

# apt update and install build dependencies
# TODO: pin package versions and check if autoconf, automake, and libtool are needed
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    g++ \
    make \
    cmake \
    git \
    wget \
    libbz2-dev \
    libcurl4-openssl-dev \
    zlib1g-dev \
    liblzma-dev \
    lbzip2 \
    autoconf \
    automake \
    libtool 
# Download and compile htslib 1.11
WORKDIR /deps
RUN wget -q https://github.com/samtools/htslib/releases/download/1.11/htslib-1.11.tar.bz2 \
    && tar -xf htslib-1.11.tar.bz2 \
    && mv htslib-1.11 htslib

WORKDIR /deps/htslib
RUN autoreconf -i && ./configure --prefix=/usr/local/ && make && make install

# Clone and compile VerifyBamID2 v2.0.1
WORKDIR /
RUN git clone --depth 1 --branch 2.0.1 https://github.com/Griffan/VerifyBamID.git

WORKDIR /VerifyBamID/build
RUN cmake .. && make && make test

# Runtime Stage using Ubuntu 24.04
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

# Install runtime dependencies
# TODO: pin package versions
RUN apt-get update && apt-get install -y --no-install-recommends \
    libgomp1 \
    libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy VerifyBamID binary and shared libs from build stage
COPY --from=build /VerifyBamID/bin/VerifyBamID /usr/local/bin/VerifyBamID2
COPY --from=build /usr/local/lib/ /usr/local/lib/

# Ensure /usr/local/lib is in the library path
ENV LD_LIBRARY_PATH=/usr/local/lib

# Make binary executable
RUN chmod +x /usr/local/bin/VerifyBamID2

# Test that the binary works
RUN VerifyBamID2 --help || (echo "VerifyBamID2 binary test failed" && exit 1)

# Default command
CMD ["VerifyBamID2"]
