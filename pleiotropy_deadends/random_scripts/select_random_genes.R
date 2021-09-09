# ---------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2020-03-04
#
# Description 
#   - Compare known genes to random set of genes
#   - and calculate a contingency table
# ---------------------------------------------------------------

# Set directory
setwd("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/data/verification_loci")

# Load packages
library(dplyr)
library(Hmisc)

# Read in gff and subset to only genes
maize_gff <- ape::read.gff("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/data/maize_v4_annotations/Zea_mays.B73_RefGen_v4.45.gff3", 
                      na.strings = c(".", "?"), GFF3 = TRUE)
maize_gff <- maize_gff[which(maize_gff$type == "gene"),]
maize_gff$geneID <- gsub("ID=gene[:]", "", maize_gff$attributes)
maize_gff$geneID <- gsub(";.*", "", maize_gff$geneID)
lala <- seq(1:10)
maize_gff <- maize_gff %>% filter(seqid %in% lala) # Only select genes on chromosomes

# Read in known genes
prior_genes <- read.csv("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/data/verification_loci/verified_genes_traits.csv",
                        header = T)
v4genes <- prior_genes$v4_gene_model

# Count number of genes within each trait category
temp <- data.frame(table(prior_genes$trait))
temp <- temp[c(11,14,4,5,9,6,17),] # Subset to traits I've mapped with

# Subtract out prior known genes from maize gff file
temp <- maize_gff[which(maize_gff$geneID %nin% v4genes), ]
temp$paper <- c("Zhang_2015", "Zhang_2015", "Zhang_2015", "Zhang_2015", "Dong_2012", "Brown_2011", "Kump_2011")

# Create a vector with the number of prior genes per chromosome
genes_to_grab <- table(prior_genes$seqid)

# Make a for loop to collect x number of genes per chromosome
random_genes <- data.frame()
for (chrom in seq(1:10)){
  genes <- maize_gff %>% filter(seqid == chrom) %>% sample_n(genes_to_grab[chrom])
  random_genes <- rbind(random_genes, genes)
}

# Get set of random genes not by chromosome
random_genes<- maize_gff %>% sample_n(nrow(prior_genes))
random_genes$trait <- prior_genes$trait

# Save random genes to file
write.csv(random_genes, "random_genes_equal_length_priors.csv", row.names = FALSE)




