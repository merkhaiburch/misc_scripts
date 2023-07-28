# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-02-02
# Updated... 2023-06-22
#
# Description:
# Conda cheetsheet
# If you change your base environment, Travis W. will find you. (unless its updating the base environment)
# ------------------------------------------------------------------------------


# ------------------------------------------------
# Environment creation and destruction
# ------------------------------------------------

# Create an environment 
conda create --name myenvironment

# Enter environment with a name
conda activate myenvironment

# Create an environment with a specific version of python
conda create --name python36-env python=3.6

# Load specific versions upon environment creation 
conda create --name basic-scipy-env ipython=7.13 matplotlib=3.1 numpy=1.10 scipy=1.4

# Create a conda environment in a specific directory
# If there are size limitations in the default conda installation path
conda create --prefix /path/to/directory/env

# Create an environment with R packages
# Preferred way to install packages, try not to do conda activate, R, install.packages()
conda create --prefix r-environment r-base r-tidyverse r-sparklyr

# Look at all environments open
conda info --envs

# Return to an environment
conda activate myenvironment

# Exit from an environment
conda deactivate

# Delete a specific environment
conda remove -n myenvironment --all

# Create a yml file (can either do manually or do these commands)
conda activate myenvironment
conda env export > environment.yml --from-history

# Create an environment from a yml file
conda env create -f environment.yml

# Update an envronmetnt (base environment shown here)
conda update --name base conda


# ------------------------------------------------
# 				  Loading Packages
# ------------------------------------------------

# Look at all installed packages in an environment (outside then inside an environment)
conda list -name myenvironment
conda list

# Install a specific package
conda install numpy

# Update a specific package
conda update numpy

# Update all packages
conda update --all

# Search for a package in the command-line (or go to the anaconda website)
conda search matplotlib

# Install a package into a pre-exisiting environment (or activate environment then install)
conda install --name myenvironment matplotlib

# Load packages from different channels (pytorch)
# Male sure paths are changed in the package list
conda create --prefix ./env --channel pytorch python=3.6 pytorch=1.5 torchvision=0.6 jupyterlab=1.0 matplotlib=3.1


# ------------------------------------------------
# 	  Wipe conda installation and start over
# ------------------------------------------------
# NOTE: using Miniconda, you'll have to specify conda-forge each time. 
# NOTE: Using miniforge is faster


# Close any versions of conda open

# Go to home directory
cd

# Show available versions of conda and where it's saving environments
conda info

# Start removing stuff (do in a screen session)
rm -rf /home/mbb262/.condarc
rm -v -rf /home/mbb262/miniconda

# For future reference (on cbsu)
rm -v -rf /home2/mbb262/.conda

# Download miniforge (linux)
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh

# close out of terminal, reopen.

# Run, answer questions and say yes. One question will ask where I want this installed, /home2/mbb262
sh miniforge3.sh

# close out of terminal, reopen 

# Check conda, update if necessary
conda info

# Install mamba, a faster solver for conda
# https://www.anaconda.com/blog/a-faster-conda-for-a-growing-community
conda install -n base conda-libmamba-solver
conda config --set solver libmamba


## For Mac
conda info
rm -rf /opt/.condarc
rm -v -rf /opt/miniconda
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh
# close out of terminal, reopen.
# Installed here
/Users/mbb262-admin/miniforge3





