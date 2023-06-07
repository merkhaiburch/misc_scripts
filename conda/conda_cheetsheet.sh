# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-02-02
# Updated... 2023-02-02
#
# Description:
# Conda cheetsheet
# ------------------------------------------------------------------------------


# ------------------------------------
# Environment creation and destruction
# ------------------------------------

# Create an environment 
conda create --name myenvironment

# Enter environment with a name
conda activate myenvironment

# Look at all environments open
conda info --envs

# Return to an environment
conda activate myenvironment

# Exit from an environment
conda deactivate

# Delete a specific environment
conda remove -n myenvironment --all


# ---------------------
# Package Stuff
# ---------------------

# Look at all installed packages in an environment (outside then inside an environment)
conda list -n myenvironment
conda list

# Install a specific package
conda install numpy

# Update a specific package
conda update numpy

# Update all packages
conda update --all


# ------------------------
# Wipe conda installation
# ------------------------

# Close any versions of conda open

# Go to home directory
cd

# Show available versions of conda
conda info

# Start removing stuff (do in a screen session)
rm -rf /home/mbb262/.condarc
rm -v -rf /home/mbb262/miniconda

# For future reference (on cbsu)
rm -v -rf /home2/mbb262/.conda

# Download miniforge
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





