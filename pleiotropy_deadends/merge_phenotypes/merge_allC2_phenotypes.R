# ---------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2020-03-03 
#
# Description 
#   - Merge C2 phenotypes
#   - Population first described in: Yang et al., 2019 PNAS
# ---------------------------------------------------------------

# C2 phenotype collector
c2_phenos <- data.frame()


# ------------------
# Yang et al., 2019
# ------------------

# Load in data
teosinte <- read.csv("~/Box Sync/Cornell_PhD/phenotype_mapping_gwas/Yang_etal_2019_PNAS/pnas.1820997116.sd01_teosinte.csv")
maize <- read.csv("~/Box Sync/Cornell_PhD/phenotype_mapping_gwas/Yang_etal_2019_PNAS/pnas.1820997116.sd02_maize.csv")

# Combine two dataframes
teosinte <- teosinte[,c(1,11:ncol(teosinte))]
maize <- maize[,c(1,12:ncol(maize))]
maize <- maize[,c(1:13,15,16,14,17:19)] # rearrange
colnames(teosinte) <- colnames(maize)
yang_2019 <- rbind(maize, teosinte)
colnames(yang_2019) <- paste0(colnames(yang_2019), "_raw_Yang2019")

# Write to CSV
write.csv(yang_2019, "all_C2_Phenos.csv", row.names = F)


