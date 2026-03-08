#!/bin/bash

rm -rf raxml-ng

# Download repo and submodules
git clone git@github.com:harmslab/raxml-ng.git
cd raxml-ng
git submodule update --init --recursive

# Build binary
mkdir build
cd build
cmake -DUSE_MPI=ON ..
make -j 4

# Copy binary to final location
cp bin/raxml-ng $CONDA_PREFIX/bin/

