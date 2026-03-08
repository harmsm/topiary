#!/bin/bash

# Download repo and submodules
git clone git@github.com:harmslab/generax.git
cd generax
git submodule update --init --recursive

# Build binary
mkdir build
cd build
cmake ..
make -j 4

# Copy binary to final location
cp bin/generax $CONDA_PREFIX/bin/

