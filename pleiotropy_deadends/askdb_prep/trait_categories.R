# ---------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2020-03-02
#
# Description 
#   - Script to make maize meeting figures
#   - Visulize trait categories
# ---------------------------------------------------------------


# Load phenotypes
all_NAM_phenos <- read.csv("~/git_projects/haplotype_v_snp_gwas/data/all_NAM_phenos.txt", sep="")
all_ames_phenos <- read.csv("~/git_projects/haplotype_v_snp_gwas/data/all_Ames_Phenos.txt", sep="")
all_goodman <- read.csv("~/git_projects/haplotype_v_snp_gwas/data/all_Assoc282_Phenos.txt", sep = "")

# Extract column names for annotation
nam <- data.frame(colnames(all_NAM_phenos))
ames <- data.frame(colnames(all_ames_phenos))
goodman <- data.frame(colnames(all_goodman))

# write.csv(, "_trait_names_annotated.csv", row.names = F)

# Set directory
setwd("~/Box Sync/Cornell_PhD/phenotype_mapping_gwas/trait_annotations/")

# Read in hand-anotated files
nam_anno <- read.csv("nam_trait_names.csv", header = TRUE)
colnames(nam_anno) <- c("trait", "category", "Population")

ames_anno <- read.csv("ames_trait_names_annotated.csv", header = T)
colnames(ames_anno) <- c("trait", "category", "Population")

goodman_anno <- read.csv("goodman_trait_names_annotated.csv", header = T)
colnames(goodman_anno) <- c("trait", "category", "Population")

c2_anno <- read.csv("c2_trait_names_annotated.csv", header = T)
colnames(c2_anno) <- c("trait", "category", "Population")

teonam_anno <- read.csv("teonam_trait_names_annotated.csv", header = T)
colnames(teonam_anno) <- c("trait", "category", "Population")

# Combine and subset
together <- rbind(nam_anno, ames_anno, goodman_anno, c2_anno, teonam_anno)
together <- together[-which(together$category == "ID"), ]
together <- together[-which(together$category == "metabolites"), ]


# Make table of counts
table_toget <- plyr::count(together, vars = c("Population", "category"))
table_toget <- table_toget[order(-table_toget$freq, table_toget$Population),]

# Change order of factor level
library(dplyr)
table_toget$name <- factor(table_toget$category, levels = as.vector(table_toget$category) %>% unique())

# Plot
library(ggplot2)
categories <- ggplot(table_toget, aes(name, freq)) +
  geom_col(aes(fill = Population)) + 
  labs(x ="Trait Category", y = "Trait x Environment Combinations") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 13), axis.title = element_text(size = 13),
        legend.position="bottom", 
        legend.text = element_text(size = 13))

ggsave(filename = "trait_categories.png", plot = categories, 
       width = 11, units = "in")


