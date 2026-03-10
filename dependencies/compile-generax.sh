#!/bin/bash

rm -rf generax

# to compile on a cluster, add a line like this one:

# module load mpi/gcc/13.1.0/openmpi/4.1.6

# This makes sure that generax compiles against the correct openmpi library. In
# a SLURM environment, you can run `module spider mpi` to find the correct
# openmpi module.

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

