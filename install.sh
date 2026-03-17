#!/bin/bash

# Check for conda
conda --version > /dev/null 2>&1
if [ $? -ne 0 ]; then
    echo "ERROR: conda not found. Please install conda or miniconda before running this script."
    exit 1
fi

# Default environment name
ENV_NAME="topiary"

# 1. Ask for environment name
echo
read -p "Enter conda environment name [$ENV_NAME]: " input_env
if [ ! -z "$input_env" ]; then
    ENV_NAME=$input_env
fi

# Check if environment already exists
conda env list | grep -q "^$ENV_NAME "
if [ $? -eq 0 ]; then
    echo
    read -p "Conda environment '$ENV_NAME' already exists. Overwrite? (y/N): " overwrite
    if [[ ! "$overwrite" =~ ^[Yy]$ ]]; then
        echo "Installation cancelled to avoid overwriting existing environment."
        exit 0
    fi
else
    overwrite="ignore"
fi

# 2. Ask for NCBI API key
echo
echo "--------------------------------------------------------------------------------"
echo "An NCBI key allows users to make more requests against the NCBI servers"
echo "when running BLAST etc. Using one is highly recommended to prevent "
echo "throttling. You may obtain one for free from the NCBI website."
echo "https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/"
echo "--------------------------------------------------------------------------------"
read -p "Enter NCBI API key (optional, press Enter to skip): " NCBI_KEY

# 3. Cluster installation prompt
echo
read -p "Are you installing on a cluster? (y/N): " is_cluster
if [[ "$is_cluster" =~ ^[Yy]$ ]]; then
    echo "--------------------------------------------------------------------------------"
    echo "To run topiary on a cluster, you need to make sure that the generax and raxml-ng"
    echo "components are compiled using your cluster's compilers and MPI. You will need"
    echo "to modify the 'dependencies/compile-generax.sh' script to do so. More details"
    echo "are at the top of the script. If you have already modified this script, please"
    echo "proceed with the installation."
    echo "--------------------------------------------------------------------------------"
    read -p "Do you want to proceed with the installation? (Y/n): " is_cluster_proceed
    if [[ "$is_cluster_proceed" =~ ^[Nn]$ ]]; then
        echo "Installation cancelled."
        exit 0
    fi
fi

# Overwrite if we have to 
if [[ "$overwrite" =~ ^[Yy]$ ]]; then
    echo "Removing existing environment '$ENV_NAME'..."
    conda deactivate 
    conda env remove -n $ENV_NAME -y
fi

# wipe local builds
rm -rf build/
for x in `find . -iname "__pycache__"`; do rm -rf $x; done

# Nuke conda environment (now handled by overwrite check above if it exists)
conda deactivate > /dev/null 2>&1

# Create new conda environment
conda env create -f environment.yml -n $ENV_NAME -y

# Set NCBI API key if provided
if [ ! -z "$NCBI_KEY" ]; then
    conda env config vars set NCBI_API_KEY=$NCBI_KEY -n $ENV_NAME
fi

# Install topiary and pip/conda dependencies
conda run -n $ENV_NAME pip install -e . -vv
conda run -n $ENV_NAME pip install coverage flake8 pytest genbadge[tests] pytest-mock sphinx pydata-sphinx-theme

# compile raxml and generax
cd dependencies
conda run -n $ENV_NAME bash compile-generax.sh
conda run -n $ENV_NAME bash compile-raxml-ng.sh

#return to start
cd ..

echo 
echo
echo "--------------------------------------------------------------------------------"
echo "Installation complete!"
echo
echo "To use topiary in the future, type 'conda activate $ENV_NAME' after opening the"
echo "terminal."
echo
if [[ "$is_cluster" =~ ^[Yy]$ ]]; then
    echo "When launching jobs on the cluster, place the lines:"
    echo
    echo "conda activate $ENV_NAME"
    echo "module load mpi/gcc/13.1.0/openmpi/4.1.6" 
    echo
    echo "at the top of your batch files. Note that the second line should match"
    echo "what you put in dependencies/compile-generax.sh"
    echo
fi
echo "--------------------------------------------------------------------------------"
echo
