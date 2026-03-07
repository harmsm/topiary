
.. include:: links.rst

.. role:: raw-html(raw)
    :format: html

.. _installation-doc:

============
Installation
============

Topiary is a python package that wraps several external software packages. This
walks through installing topiary and all of the external software packages it
wraps. Note that topiary requires a linux or macOS machine.

The basic installation instructions for topiary are:

1. Download and install `miniconda <miniconda-link_>`_
2. Install topiary using conda

1. Install miniconda
====================

To prevent interference with other packages, we recommend installing topiary in
its own conda environment. If you do not have conda installed, download and
install `miniconda <miniconda-link_>`_ before proceeding. Installation instructions
are available on the `linked miniconda page <miniconda-link_>`_. 

2. Install topiary
==================

Once conda is installed, open a standard terminal. Copy the following commands
into the prompt to run them. 

.. code-block:: shell-session

  conda install git
  git clone https://github.com/harmslab/topiary
  cd topiary
  conda config --set channel_priority strict
  conda env create -f environment.yml
  conda activate topiary
  python -m pip install . -vv

At this point, you should have topiary and its supporting software installed. 

2.1. Installing RAxML-NG and GeneRax
====================================

Unfortunately, the versions of RAxML-NG and GeneRax on conda are currently
broken. To get around this, we have to download and compile these software
packages manually. In a shell, run:

.. code-block:: shell-session

  git clone https://github.com/harmsm/raxml-ng.git
  cd raxml-ng
  mkdir build
  cd build
  cmake ..
  make -j 4
  cp build/bin/raxml-ng $CONDA_PREFIX/bin/

.. code-block:: shell-session

  git clone https://github.com/harmsm/generax.git
  cd generax
  mkdir build
  cd build
  cmake ..
  make -j 4
  cp build/bin/generax $CONDA_PREFIX/bin/


3. Check supporting software packages
=======================================

You can check which software packages are visible to topiary by:

.. code-block:: shell-session

  topiary-check-installed

The output should look something like this:

.. image:: _static/img/installation/topiary-check-installed_450x483.png
  :align: center
  :alt: topiary-check-installed terminal output

:raw-html:`<br />`
If some of the packages are not installed (:code:`passes: N`), proceed to the
sections below.

----------------------------
NCBI API Key
----------------------------

If you wish to use an `NCBI API Key <https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/>`_,
set the environment variable :code:`NCBI_API_KEY` to point to your key. For 
example: 

.. code-block:: shell-session

  export NCBI_API_KEY='abcdef012'

Topiary will recognize this key and increase the number of allowed requests per
second to NCBI servers. 

If any of the packages are not installed, you can try to manually install then with:

.. code-block:: shell-session

  conda install -c conda-forge -c bioconda mpi4py openmpi "muscle>=5.0" "raxml-ng>=1.1" "generax>=2.0" "blast>=2.2"

You can then check to make sure everything installed correctly by running:

.. code-block:: shell-session

  topiary-check-installed

If any of these packages were not installed by conda--or you wish to install
them yourself--you can install them manually using the following links:

+ `NCBI blast+ >= 2.2 <blast-download_>`_. (This will install both the blastp and
  makeblastdb programs.)
+ `muscle >= 5.0 <muscle-download_>`_.
+ `GeneRax >= 2.1.3 <generax-download_>`_.
+ `RAxML-NG >= 1.2.2 <raxml-ng-download_>`_.

After installation, you'll need to make sure the directories containing these
binaries are in your :code:`$PATH` directory. (See `here <nix-path_>`_ for
instructions).

------------------
Required libraries
------------------

+ Core scientific python libraries:

  + `Python >= 3.8 <python-link_>`_
  + `numpy <numpy-link_>`_
  + `pandas <pandas-link_>`_
  + `matplotlib <matplotlib-link_>`_

+ Tree manipulation/drawing:

  + `ete3 <ete3-download_>`_
  + `toytree <toytree-download_>`_
  + `dendropy <dendropy-download_>`_

+ Packages used for tree/ancestor inferences:

  + `NCBI BLAST+ <blast-download_>`_
  + `muscle >= 5.0 <muscle-download_>`_
  + `GeneRax >= 2.0 <generax-download_>`_
  + `RAxML-NG >= 1.1 <raxml-ng-download_>`_
  + `python-opentree <opentree-link_>`_
