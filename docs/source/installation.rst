
.. include:: links.rst

.. role:: raw-html(raw)
    :format: html

.. _installation-doc:

============
Installation
============

Topiary is a python package that wraps several external software packages. This
walks through installing topiary and all of the external software packages it
wraps. 

.. IMPORTANT::
   Topiary requires a linux or macOS machine. **Windows is no longer supported.**

Basic Installation
==================

The easiest way to install topiary and its dependencies is to use the provided
installation script. 

1. Download and install `miniconda <miniconda-link_>`_
2. Open a terminal and run the following commands:

.. code-block:: shell-session

  conda install git
  git clone https://github.com/harmslab/topiary
  cd topiary
  bash install.sh

The `install.sh` script will:

* Prompt you for a conda environment name (default is `topiary`).
* Ask for an optional `NCBI API Key <https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/>`_. Providing one is highly recommended to avoid rate-limiting during BLAST searches.
* Prompt you about cluster installation (see below).
* Create the conda environment, install all python dependencies, and compile the necessary external binaries (RAxML-NG and GeneRax).

Cluster Installation
====================

If you are installing topiary on a high-performance computing (HPC) cluster, you
generally need to ensure that the external software (GeneRax) is
compiled using the cluster's specific MPI and C compilers. 

To do this:

1. Identify the correct MPI and compiler modules on your cluster (e.g., using `module spider mpi`). The details of this step depend on your cluster architecture. If you do not know how to do this, please contact your cluster administrator.
2. Open `dependencies/compile-generax.sh` in a text editor.
3. Add the appropriate `module load` command near the top of the script. For example:

   .. code-block:: bash

      module load mpi/gcc/13.1.0/openmpi/4.1.6

4. Save the file and run `bash install.sh`. When prompted if you are on a cluster, answer `y`. 

This ensures that GeneRax and RAxML-NG are linked against the cluster's high-performance interconnects, allowing them to run efficiently across multiple nodes.

Checking Installation
=====================

After the installation script finishes, you can verify that everything is 
installed correctly by running:

.. code-block:: shell-session

  conda activate topiary
  topiary-check-installed

(Replace `topiary` with your custom environment name if you chose one).

The output should show :code:`passes: Y` for all required packages:

.. image:: _static/img/installation/topiary-check-installed_450x483.png
  :align: center
  :alt: topiary-check-installed terminal output

----------------------------
NCBI API Key
----------------------------

If you did not provide an NCBI API key during the initial installation, you can 
add one later by running:

.. code-block:: shell-session

  conda env config vars set NCBI_API_KEY='your_key_here' -n topiary

This ensures the key is automatically exported whenever the topiary environment 
is activated.

Required Software
=================

While `install.sh` handles these automatically, topiary relies on the following
software. Note that we use custom versions of these (specifically GeneRax and
RAxML-NG) so we cannot gaurantee topiary will function using other versions of 
these pieces of code. Our custom versions are included in the `dependencies`
directory.

+ `NCBI blast+ >= 2.2 <blast-download_>`_
+ `muscle >= 5.0 <muscle-download_>`_
+ `GeneRax >= 2.1.3 <generax-download_>`_
+ `RAxML-NG >= 1.2.2 <raxml-ng-download_>`_
+ `Python >= 3.11 <python-link_>`_

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
  + `GeneRax >= 2.1.3 <generax-download_>`_
  + `RAxML-NG >= 1.2.2 <raxml-ng-download_>`_
  + `python-opentree <opentree-link_>`_
