# ---------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2019-12-06
#
# Description 
# 	- Convert gene IDs between any maize version
# ---------------------------------------------------------------

# v3 to v4 translation
v3_v4_xref <- read.delim("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/data/maize_v4_annotations/v3_v4_xref.txt")


# -----------------------------------------
# Format Dong 2012 gene list from Brandon
# One time use only
# -----------------------------------------

# Load package
# BiocManager::install("genomes")
library(ape)

# Load in Dong 2012 genes formatted from
# https://bitbucket.org/bucklerlab/askdb_cis_trans/src/master/data/flowering_genes_dong.csv
dong <- read.csv("flowering_genes_dong.csv", header = TRUE)

# Read gff file
maize_gff <- ape::read.gff("~/Box\ Sync/Cornell_PhD/labProjects/hap_gwa/ames_tests/Sep25_localPCs/Zea_mays.B73_RefGen_v4.45.gff3" , na.strings = c(".", "?"), GFF3 = TRUE)

# Do some formatting to isolate gene names
maize_gff <- maize_gff[which(maize_gff$type == "gene"),]
maize_gff$gene <- gsub("ID=gene:", "", maize_gff$attributes)
maize_gff$gene <- gsub(";.*", "", maize_gff$gene)

# Merge dong genes with gff annotations
dong_coordinates <- merge(x = dong, y = maize_gff,
                          by.x = "v4_assoc_gene_model", by.y = "gene",
                          all.x = TRUE) 
dong_coordinates <- dong_coordinates[,c(2:8,1,9, 12:13,17)]

# Write csv
write.csv(dong_coordinates, "dong_genes_with_coordinates.csv", row.names = F, quote = T)


# -------------------------------------
# Format Zhang 2015
# priori genes for C and N metabolism
# Taken from Supplemental Figures 1-13 and Supplemental Tables 1 and 8
# https://doi.org/10.1104/pp.15.00025
# -------------------------------------

# Read in file
zhang <- read.csv("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/data/Zhang2014_metabolism_priors.csv")

# Read in gene model association file (conversion file) and do formatting
gene_model_xref_v3 <- read.delim("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/data/maize_v4_annotations/gene_model_xref_v3.txt")

# Merge my genes with annotations (v1 --> v4)
zhang_allModels <- merge(x = zhang, y = gene_model_xref_v3,
                         by.x = "gene", by.y = "v1_gene_model")

# Subset out columns
zhang_allModels <- zhang_allModels[,c(14,3,2,5,4,12)]

# Merge my genes with gff annotations from v4
zhang_gene_coords <- merge(x = zhang_allModels, y = maize_gff,
                           by.x = "v4_gene_model", by.y = "gene")

# Merge with dong flowering time loci
# Get data
dong_genes_with_coordinates <- read.csv("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/ames_tests/Sep25_localPCs/dong_genes_with_coordinates.csv", header = TRUE)

# Fill out missing gaps between datasets
dong_genes_with_coordinates$Paper <- rep("dong_2012", nrow(dong_genes_with_coordinates))
dong_genes_with_coordinates$trait <- rep("flowering_time", nrow(dong_genes_with_coordinates))
zhang_gene_coords$pathway <- rep(NA, nrow(zhang_gene_coords))
# zhang_gene_coords$zm_gene <- rep(NA, nrow(zhang_gene_coords))
zhang_gene_coords$at_gene <- rep(NA, nrow(zhang_gene_coords))
zhang_gene_coords$tissues_expressed <- rep(NA, nrow(zhang_gene_coords))
zhang_gene_coords$genbank_accession <- rep(NA, nrow(zhang_gene_coords))

# Rearrange my data
zhang_gene_coords <- zhang_gene_coords[,c(1,4,7,6,10,11,3,5,17,18,19,15)]
dong_genes_with_coordinates <- dong_genes_with_coordinates[,c("v4_assoc_gene_model","v3_assoc_gene_model","seqid",
                                                              "zm_gene","start","end","trait","Paper",
                                                              "at_gene","tissues_expressed","genbank_accession",
                                                              "attributes")]
# Make colnames consistent
colnames(dong_genes_with_coordinates) <- colnames(zhang_gene_coords)

# Combine
prior_genes <- rbind(dong_genes_with_coordinates, zhang_gene_coords)

# Write file
write.csv(prior_genes, "verified_genes_traits.csv", quote = F, row.names = F)

