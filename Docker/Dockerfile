# Define base image Ubuntu 20.04 as Build Stage
FROM ubuntu:20.04 AS build

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive
ENV LD_LIBRARY_PATH=/usr/local/lib

# apt update and install build dependencies

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    g++=4:9.3.0-1ubuntu2 \
    make \
    cmake \
    git \
    wget \
    libbz2-dev=1.0.8-2 \
    libcurl4-openssl-dev \
    zlib1g-dev \
    liblzma-dev \
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

RUN apt-get update && apt-get install -y --no-install-recommends \
    libgomp1=14.2.0-4ubuntu2~24.04 \
    libcurl4-openssl-dev=8.5.0-2ubuntu10.6 \
    && rm -rf /var/lib/apt/lists/*

# Copy VerifyBamID binary and shared libs from build stage
COPY --from=build /VerifyBamID/bin/VerifyBamID /usr/local/bin/VerifyBamID2
COPY --from=build /usr/local/lib/ /usr/local/lib/

# Ensure /usr/local/lib is in the library path
ENV LD_LIBRARY_PATH=/usr/local/lib

# Make binary executable
RUN chmod +x /usr/local/bin/VerifyBamID2

# Test that the binary works
RUN VerifyBamID2 --help 2>&1 | grep -q "^ Version:" || (echo "VerifyBamID2 binary test failed" && exit 1)

ENTRYPOINT ["/usr/local/bin/VerifyBamID2"]
