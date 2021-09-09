# ---------------------------------------
# Merritt Burch
# mbb262@cornell.edu
# 2019-07-19
# Script to merge data together from various
# QTL/gwas studies and appending the author
# name to the trait so they don't get lost
# This script merges all IBM phenotypes together
# ---------------------------------------

ibm_phenos <- data.frame()

# Full IBM names
ibm <- read.table("~/Box Sync/Cornell_PhD/phenotype_mapping_gwas/ibm_taxa.txt",
           header = TRUE)
