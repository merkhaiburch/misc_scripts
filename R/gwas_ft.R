# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-05-14 
# Updated... 2022-05-14

# Description 
# Check mapping of traits in the 282 for QC reasons
# ---------------------------------------------------------------

# Load in packages
library(data.table)
library(dplyr)
library(ggplot2)


# Try looking at buckler phenotypes
flowering_phenos <- data.table::fread("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/data/phenos/2021_09_final_update/goodman/02_goodman_phenotypes_Buckler_2009.txt")
colnames(flowering_phenos)[1] <- "Trait"
cor(flowering_phenos[,2], flowering_phenos[,3])

flowering_phenos_2 <- data.table::fread("~/git_projects/pleiotropy/data/all_Assoc282_Phenos.csv")
flowering_phenos_2 <- flowering_phenos_2 %>% select(Geno_Code, hapmap_names, DTA_BLUP_Buckler_2009_goodman, DTS_BLUP_Buckler_2009_goodman, ASI_BLUP_BLUP_Buckler_2009_goodman)
temp <- merge(flowering_phenos, flowering_phenos_2, by.x = "Trait", by.y = "hapmap_names")

# Load in original file
og_flowering <- data.table::fread("~/Box Sync/Cornell_PhD/phenotype_mapping_gwas/buckler_2009_science_wholeplant_floweringtime/NAMSum0607FloweringTraitBLUPsAcross8Envs.csv")
colnames(og_flowering)[1] <- "Trait"

# Merge two files together
temp2 <- merge(temp, og_flowering, by.x = "Geno_Code", by.y = "Trait")

# The phenotypes are the same, the problem isn't here!!
cor(temp2[,3], temp2[,13])
cor(temp2[,4], temp2[,14])
cor(temp2[,15], temp2[,15])


# -------------------------------------------------------------------------------------------

# Load in flowering time results from buckler 2009 NAM 
buckler_files <- list.files("~/Box Sync/Cornell_PhD/labProjects/debugging/",
                            pattern = "*filtered.csv", full.names = TRUE)
buckler2009 <- rbindlist(lapply(buckler_files, data.table::fread))

# Subset into dta, dts, and asi results
dta <- buckler2009 %>% filter(Trait == "DTA_BLUP_Buckler_2009_goodman")
dts <- buckler2009 %>% filter(Trait == "DTS_BLUP_Buckler_2009_goodman")
asi <- buckler2009 %>% filter(Trait == "ASI_BLUP_BLUP_Buckler_2009_goodman")



