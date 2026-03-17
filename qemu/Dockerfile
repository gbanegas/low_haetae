FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

# Install build tools, ARM cross toolchain and QEMU
RUN apt-get update && apt-get install -y \
    build-essential \
    make \
    git \
    gcc-arm-none-eabi \
    binutils-arm-none-eabi \
    libnewlib-arm-none-eabi \
    qemu-system-arm \
    gdb-multiarch \
 && rm -rf /var/lib/apt/lists/*

# Work directory inside the container
WORKDIR /workspace

# Default: drop into a shell
CMD ["bash"]
