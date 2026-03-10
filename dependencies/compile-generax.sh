#!/bin/bash

rm -rf generax

# to compile on a cluster, add a line like this one. In a SLURM environment, 
# run module spider mpi to find the correct openmpi/gcc binaries.
#module load mpi/gcc/13.1.0/openmpi/4.1.6

# Download repo and submodules
git clone git@github.com:harmslab/generax.git
cd generax
git submodule update --init --recursive

# Build binary
mkdir build
cd build
cmake -DUSE_MPI=ON ..
make -j 4

# Copy binary to final location
cp bin/generax $CONDA_PREFIX/bin/

