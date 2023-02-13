# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-02-02
# Updated... 2023-02-02
#
# Description:
# Conda cheetsheet
# ------------------------------------------------------------------------------

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

# Look at all installed packages in an environment (outside then inside an environment)
conda list -n myenvironment
conda list