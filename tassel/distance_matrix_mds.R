# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-04-19 
# Updated... 2022-04-19

# Description 
# Calculate distance matrix, calculate MDS on distance matrix
# ---------------------------------------------------------------

# set memory
options(java.parameters = c("-Xmx50g"))

# Load packages
library(dplyr)
library(rTASSEL)

# Load in test genotype data
genoPathVCF <-  rTASSEL::readGenotypeTableFromPath("~/Desktop/TutorialData/mdp_genotype.hmp.txt",
                                           keepDepth = FALSE)
print(genoPathVCF)