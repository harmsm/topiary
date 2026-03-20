#!/bin/bash

# NOTE: this script should generally be run by "install.sh", not run by itself.

# -----------------------------------------------------------------------------
# COMPILING ON A CLUSTER
# 
# To run generax on a cluster, you need to make sure that it is compiled
# using your cluster's compilers and MPI. To do this, you can add a line like
# this one to the top of this script:

# module load mpi/gcc/13.1.0/openmpi/4.1.6

# This makes sure that generax compiles against the correct openmpi library, 
# which allows it to run across multiple cores/nodes on a cluster. The exact
# command depends on your cluster. If your cluster uses SLURM, you can usually
# find the correct command by running `module spider mpi`. If you are unsure
# which command to use, please contact your cluster administrator.

# -----------------------------------------------------------------------------
# No user modifications should be required below this line

# Check for conda environment
if [ -z "$CONDA_PREFIX" ]; then
    echo "ERROR: \$CONDA_PREFIX is not set. This script should be run via "
    echo "'conda run -n <env_name> bash compile-generax.sh'"
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
    if [ -f "$CONDA_PREFIX/bin/generax" ]; then
        echo "Using existing generax binary in \$CONDA_PREFIX/bin/ (skip compile)"
        exit 0
    fi
fi

# Clean up environment
rm -rf generax

# Download repo and submodules
git clone https://github.com/harmslab/generax.git
cd generax
git submodule update --init --recursive

# Build binary
mkdir build
cd build
cmake -DUSE_MPI=ON -DCORAX_BUILD_PORTABLE=ON ..
make -j 4

# Copy binary to final location
cp bin/generax $CONDA_PREFIX/bin/

