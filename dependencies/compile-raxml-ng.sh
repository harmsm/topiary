#!/bin/bash

# This script should generally be run by "install.sh", not run by itself. It 
# should not require any user edits. 

# Check for conda environment
if [ -z "$CONDA_PREFIX" ]; then
    echo "ERROR: \$CONDA_PREFIX is not set. This script should be run via "
    echo "'conda run -n <env_name> bash compile-raxml-ng.sh'"
    exit 1
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
cmake ..
make -j 4

# Copy binary to final location
cp bin/raxml-ng $CONDA_PREFIX/bin/

