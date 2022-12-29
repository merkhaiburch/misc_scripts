# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-09-22 
# Updated... 2022-09-22
#
# Description:
# Start a kotlin kernel in a conda environment
# ------------------------------------------------------------------------------

# install the kotlin kernel and jupyter lab on cbsu machine
conda create -c conda-forge -c jetbrains -n evanEnv jupyterlab kotlin-jupyter-kernel

# Activate my conda environment
conda activate merrittEnv

# Start Jupyter notebook server, go to the link it provides in a vnc session
jupyter lab

# To deactivate conda session
conda deactivate