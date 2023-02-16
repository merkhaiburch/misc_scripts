# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-02-01
# Updated... 2023-02-01
#
# Description:
# investigate expressed genes in example b73, mo17, and hybrid b73xmo17 alignment
# ------------------------------------------------------------------------------

# Load packages
library(data.table)
library(dplyr)

# Load in b73 counts
b73 <- data.table::fread("~/Downloads/b73_counts.txt", header = TRUE)
colnames(b73) <- gsub("_Zm.*", "", colnames(b73))

# load in mo17 counts
mo17 <- data.table::fread("~/Downloads/mo17_counts.txt", header = TRUE)
colnames(mo17) <- gsub("_Zm.*", "", colnames(mo17))

# Load in hybrid counts
hybrid <- data.table::fread("~/Downloads/hybrid_b73_mo17_counts.txt", header = TRUE)
# colnames(hybrid) <- gsub("_Zm.*", "", colnames(hybrid))

# Look at dim
b73[1,1:5]
mo17[1,1:5]
hybrid[1,1:5]

#Format dfs
b73 <- t(b73) %>% as.data.frame()
colnames(b73) <- b73[1,]
b73 <- b73[-1,]
b73[,1] <- as.numeric(b73[,1])
b73[,2] <- as.numeric(b73[,2])
str(b73)

mo17 <- t(mo17) %>% as.data.frame()
colnames(mo17) <- mo17[1,]
mo17 <- mo17[-1,]
mo17[,1] <- as.numeric(mo17[,1])
mo17[,2] <- as.numeric(mo17[,2])
str(mo17)

hybrid <- t(hybrid) %>% as.data.frame()
colnames(hybrid) <- hybrid[1,]
hybrid <- hybrid %>% filter(hybridB73Mo17_to_jointB73Mo17_minimap2 != "hybridB73Mo17_to_jointB73Mo17_minimap2")
hybrid[,1] <- as.numeric(hybrid[,1])
str(hybrid)
hybrid$taxa <- gsub(".*_", "", rownames(hybrid))
hybrid$gene <- gsub("_.*", "", rownames(hybrid))

# Count how many genes are expressed
b73 %>% 
  filter(b73_to_b73_minimap2 > 0) %>% 
  summarise(n = n()) # 24399

b73 %>% 
  filter(hybridB73Mo17_to_b73_minimap2 > 0) %>% 
  summarise(n = n()) # 23122

mo17 %>% 
  filter(hybridB73Mo17_to_mo17_minimap2 > 0) %>% 
  summarise(n = n()) # 21429

mo17 %>% 
  filter(mo17_to_mo17_minimap2 > 0) %>% 
  summarise(n = n()) # 21346

hybrid %>% 
  group_by(taxa) %>% 
  filter(hybridB73Mo17_to_jointB73Mo17_minimap2 > 0) %>% 
  summarise(n = n()) # b73 = 14937, mo17 = 12818


# Check correlations between individual genes in hybrids, are they similar?






