
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "This script should be run as 'source fresh-start.sh', not 'bash fresh-start.sh'."
    exit
fi

# wipe local builds
rm -rf build/
for x in `find . -iname "__pycache__"`; do rm -rf $x; done

# Nuke conda environment
conda deactivate
conda env remove -n topiary -y

# Create new conda environment
conda env create -f environment.yml -n topiary -y

# Install topiary and pip/conda dependencies
conda run -n topiary pip install -e . -vv
conda run -n topiary pip install coverage flake8 pytest genbadge[tests] pytest-mock
conda activate topiary

# compile raxml and generax
cd dependencies
bash compile-generax.sh
bash compile-raxml-ng.sh

#return to start
cd ..
