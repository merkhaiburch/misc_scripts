# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-07-27 
# Updated... 2022-07-27

# Description 
# Gather unique taxa names from NAM to subset Hapmap files with
# ---------------------------------------------------------------

# Load packages
library(magrittr)
library(dplyr)

# Load in metadata file
taxa <- read.csv("~/Box Sync/Cornell_PhD/labProjects/nam_hybrid_te/data/sample_metadata.csv")

# Isolate unique taxa
temp <- taxa %>% filter(type == "inbred") %$% unique(cultivar)

# Subset the inbred names based on matches --> did this one by hand
inbred <- read.delim("~/Box Sync/Cornell_PhD/labProjects/nam_hybrid_te/data/te_ase_nam_inbreds_hapmapids.txt", header = F)

# Create a file with B73 x NAM inbreds, PHZ51 x NAM inbreds, and PHB47 x NAM inbreds
b73 <- inbred
b73$p1 <- rep("282set_B73", nrow(b73))
b73 <- b73[-1,] # remove first line which is just b73

phz51 <- inbred
phz51$p1 <- rep("PHZ51", nrow(phz51))
phz51 <- phz51[-28,]

b47 <- inbred
b47$p1 <- rep("B47", nrow(b47))
b47 <- b47[-29,]

# combine
hybrids <- rbind(b73, phz51, b47)

# export for use in tassel's Create Hybrid Genotypes
write.table(hybrids, "~/Box Sync/Cornell_PhD/labProjects/nam_hybrid_te/data/te_ase_hybrids_for_tassel.txt", 
          quote = F, row.names = F, sep = "\t", col.names = F)

