#!/bin/bash

# This script should generally be run by "install.sh", not run by itself. It 
# should not require any user edits. 

# Check for conda environment
if [ -z "$CONDA_PREFIX" ]; then
    echo "ERROR: \$CONDA_PREFIX is not set. This script should be run via "
    echo "'conda run -n <env_name> bash compile-raxml-ng.sh'"
    exit 1
fi

KEEP_EXISTING=0
while [[ $# -gt 0 ]]; do
  case $1 in
    --keep-existing)
      KEEP_EXISTING=1
      shift ;;
    *)
      shift ;;
  esac
done

# Check if binary already exists in final location
if [ $KEEP_EXISTING -eq 1 ]; then
    if [ -f "$CONDA_PREFIX/bin/raxml-ng" ]; then
        echo "Using existing raxml-ng binary in \$CONDA_PREFIX/bin/ (skip compile)"
        exit 0
    fi
fi

# Clean up environment
rm -rf raxml-ng

# Download repo and submodules
git clone https://github.com/harmslab/raxml-ng.git
cd raxml-ng
git submodule update --init --recursive

# Build binary
mkdir build
cd build
cmake -DTERRAPHAST_ARCH_NATIVE=OFF -DENABLE_RAXML_SIMD=OFF -DENABLE_PLLMOD_SIMD=OFF -DENABLE_SSE=OFF -DENABLE_AVX=OFF -DENABLE_AVX2=OFF ..
make -j 4

# Copy binary to final location
cp bin/raxml-ng $CONDA_PREFIX/bin/

