# ---------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2020-02-19
#
# Description 
#   - Parameter tuning for GWAS models
#   - Trying many models
# 
# ---------------------------------------------------------------

# Load in source scripts
source('~/git_projects/haplotype_v_snp_gwas/src/R/pca_testing/ames2any_matrix.R')
source('~/git_projects/haplotype_v_snp_gwas/src/R/pca_testing/local_window_funs.R')
source("~/git_projects/haplotype_v_snp_gwas/src/R/pca_testing/count_gwas_overlaps.R")
source('~/git_projects/haplotype_v_snp_gwas/src/R/random_scripts/manhattan_plot.R')

# Load packages
library(gdsfmt)
library(SNPRelate)
library(ggplot2)
library(ape)
library(dplyr)
library(broom)
library(magrittr)
library(patchwork)


# ------------------
#   Load test data
# ------------------

# Load in Beagle imputed SNPs that were filtered with bcftools
ames.vcf <- "~/Box Sync/Cornell_PhD/labProjects/hap_gwa/data/genos/hapmap_v321_snps/ames_for_pca_merged.vcf"
ames_vcf <- snpgdsVCF2GDS(ames.vcf, "ames.gds", method="copy.num.of.ref", snpfirstdim=TRUE, ignore.chr.prefix = "S")
snpgdsSummary("ames.gds")
genofile_ames <- snpgdsOpen(ames_vcf)

# Load in NAM SNPs that are shared (intersect) with the Ames SNPs
# nam.vcf <- "/workdir/mbb262/genotypes/nam/imputed/combined_ibm_nam_all_imputed_bcftools/filtered_for_pca"
nam.vcf <- "~/Box Sync/Cornell_PhD/labProjects/hap_gwa/data/genos/hapmap_v321_snps/nam_by_ames_sites_for_pca_allChroms.vcf"
nam_vcf <- snpgdsVCF2GDS(nam.vcf, "nam.gds", method="copy.num.of.ref", snpfirstdim=TRUE, ignore.chr.prefix = "S")
snpgdsSummary("nam.gds")
genofile_nam <- snpgdsOpen(nam_vcf)

# Get positions of all snps within snps
temp <- snpgdsSNPList(genofile_nam)
chr10_temp <- temp %>% filter(chromosome == 10, position >= 141561862, position <=141564767)
plot(chr10_temp$snp.id, chr10_temp$position)

plot(temp$snp.id, temp$position)
# ------------------------------------------
#     Calculate ames2nam global PCs
#     uses ames2any_matrix functions
# ------------------------------------------

# Combine genofile objects
X_allNAM <- combine_genofile_matrix(genofile_ames = genofile_ames, genofile_any = genofile_nam)

# ------
# 3 gPCs
# ------ 

# Calculate Ames global PCs, get loadings
global_ames_BV_3gPCs <- get_ames_BV(X_ames = genofile_ames, Q = NULL, num_PCs = 3)

# Transfer loadings and coefficients over to all X (Ames, NAM) to get adjusted global PCs
# Get Q matrix (population structure)
ames2nam_gPCs_3gPCs_allNAM <- ames2any(X_all = X_allNAM$X_all,
                                       Q = NULL, B = global_ames_BV_3gPCs$B, V = global_ames_BV_3gPCs$V)
Q_all_3gPCs_allNAM <- cbind(1, ames2nam_gPCs_3gPCs_allNAM)

# ------
# 26 gPCs
# ------ 

# Calculate Ames global PCs, get loadings
global_ames_BV_26gPCs <- get_ames_BV(X_ames = genofile_ames, Q = NULL, num_PCs = 26)
ames2nam_gPCs_26gPCs_allNAM <- ames2any(X_all = X_allNAM$X_all,
                                        Q = NULL, B = global_ames_BV_26gPCs$B, V = global_ames_BV_26gPCs$V)
Q_all_26gPCs_allNAM <- cbind(1, ames2nam_gPCs_26gPCs_allNAM)


# ------------------------------------------
#   Calculate local PCs in defined windows
#           formatting step
# ------------------------------------------

# On desktop this file is typically in
# "~/Box Sync/Cornell_PhD/labProjects/hap_gwa/data/maize_v4_annotations/Zea_mays.B73_RefGen_v4.45.gff3"

# Get all maize v4 genes from gff file 
maize_gff <- ape::read.gff("/workdir/mbb262/data/Zea_mays.B73_RefGen_v4.45.gff3", 
                      na.strings = c(".", "?"), GFF3 = TRUE)
maize_gff <- ape::read.gff("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/data/maize_v4_annotations/Zea_mays.B73_RefGen_v4.45.gff3", 
                      na.strings = c(".", "?"), GFF3 = TRUE)

# Do some formatting to isolate gene names
maize_gff <- maize_gff[which(maize_gff$type == "gene"),]
maize_gff$gene <- gsub("ID=gene:", "", maize_gff$attributes)
maize_gff$gene <- gsub(";.*", "", maize_gff$gene)

# Rearrange, get rid of columns I don't need
genes <- maize_gff[, c(10,1,4,5)]
genes <- genes %>% 
  dplyr::filter(seqid %in% 1:10) %>% 
  dplyr::mutate(seqid = factor(.$seqid, levels = 1:10))

# Use the function to cut genome into 150 gene windows
# 177 genes = 226 windows; 178 = 225; 179=221; 180=220.33; 181 = 220!
window_file <- getWindowCoords(362, genes)

# Format file
window_file <- as.data.frame(window_file)
colnames(window_file) <- c("window", "chrom", "start", "stop")
window_file$chrom <- as.numeric(as.character(window_file$chrom))
window_file$start <- as.numeric(as.character(window_file$start))
window_file$stop <- as.numeric(as.character(window_file$stop))


# -------------------------------------------
# Local PC function with gene number windows
# -------------------------------------------

# Use on variable gene number windows with ames to nam PCs
# 3 gPCs
geneWindow_ames2nam_local_3allNAM <- calc_local_pcs_gene_window(X_all = X_allNAM,
                                                                genofileID_ames = genofile_ames,
                                                                Q_all = Q_all_3gPCs_allNAM,
                                                                numLocalPcs = 3,
                                                                windowFile = window_file)

# 26 gPCs
geneWindow_ames2nam_local_26allNAM <- calc_local_pcs_gene_window(X_all = X_allNAM,
                                                                 genofileID_ames = genofile_ames,
                                                                 Q_all = Q_all_26gPCs_allNAM,
                                                                 numLocalPcs = 3,
                                                                 windowFile = window_file)


# -----------------------------
# Export files in tassel format
# -----------------------------

# Load in NAM data
# Filter down phenotypes to only 
# fumarate,malate,chlorophyll A,chlorophyll b, cobdiameter, southern leaf blight, flowering time
# File is typically in
# "~/Box Sync/Cornell_PhD/labProjects/hackathons/2019-10-21_hackathon/all_NAM_phenos.txt"
# all_NAM_phenos <- read.csv("/workdir/mbb262/data/all_NAM_phenos.txt", sep="")
# all_NAM_phenos <- read.csv("~/Box Sync/Cornell_PhD/labProjects/hackathons/2019-10-21_hackathon/all_NAM_phenos.txt", sep = "")
all_NAM_phenos <- read.table("~/git_projects/haplotype_v_snp_gwas/data/all_NAM_phenos.txt", sep="", header = TRUE)
all_NAM_phenos <- all_NAM_phenos[,c(2,4,40,39,37,17,8:9)]

# "Fix" IBM names
all_NAM_phenos$IBM_compat <- as.character(all_NAM_phenos$Geno_Code)
all_NAM_phenos[c(3203:3402),9] <- as.character(all_NAM_phenos[c(3203:3402),2])
all_NAM_phenos$IBM_compat <- gsub("MO", "M0", all_NAM_phenos$IBM_compat)

# Function to format data
tasselize <- function(main_matrix, global_matrix, window_matrix, filename){
  # Combine objects
  together <- cbind(data.frame(rownames(main_matrix)), data.frame(global_matrix), data.frame(window_matrix))
  
  # Format taxa IDs from PCs for merging
  # together$copied_taxa <- gsub(":[0-9]{9}", "", together[,1])
  # together$copied_taxa <- gsub("_", "", together$copied_taxa)
  together$copied_taxa <- together[,1]
  
  # Merge dataframes together
  together <- merge(x = all_NAM_phenos, y = together, sort = F,
                    by.x = "IBM_compat", by.y = "copied_taxa")
  
  # Sort columns (3 or 26 gPCs)
  together <- together[,c(13,4:12,14:ncol(together))]
  
  # Change name of first column from taxa ID to trait
  colnames(together)[1] <- "<Trait>"
  
  # Remove missing and duplicated values
  together <- na.omit(together)
  together <- together[!duplicated(together),]
  
  # To export dataset as is to work on in R
  return(together)
  
  # Subset out only the global PCs
  global <- together[,c(1,11:13)]
  global[,1] <- as.character(global[,1])
  header_global <- c("<Covariate>", rep("", ncol(global)-1))
  global <- rbind(header_global, colnames(global), global)
  
  # Subset out phenotypes and format
  phenos <- together[,1:10]
  phenos[,1] <- as.character(phenos[,1])
  header_phenos <- c("<Phenotype>", rep("", ncol(phenos)-1))
  phenos <- rbind(header_phenos, colnames(phenos), phenos)
  
  # Make column names the first row
  together[,1] <- as.character(together[,1])
  header <- c("<Covariate>", rep("", ncol(together)-1))
  together <- rbind(header, colnames(together), together)
  
  # Only get PCs
  together <- together[,c(1,11:ncol(together))]
  
  # Export
  write.table(together, file = filename, row.names = FALSE, 
              col.names = FALSE, sep = "\t", quote = FALSE)
  
  # Export phenotypes
  write.table(phenos, file = "wallace_phenos_2020_03_10.txt", row.names = FALSE,
              col.names = FALSE, sep = "\t", quote = FALSE)
  
  # Export global PCs only
  # write.table(global, file = "ames2nam_3gPCs_allNAM_gPCs_only.txt",
  #             row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  
}


# Gene window options
tasselize(X_allNAM$X_all, ames2nam_gPCs_3gPCs_allNAM, geneWindow_ames2nam_local_3allNAM, 
          filename = "ames2nam_3gPCs_allNAM_geneWindow.txt")
tasselize(X_allNAM$X_all, ames2nam_gPCs_26gPCs_allNAM, geneWindow_ames2nam_local_26allNAM, 
          filename = "ames2nam_26gPCs_allNAM_geneWindow.txt")

temp <- tasselize(X_allNAM$X_all, ames2nam_gPCs_3gPCs_allNAM, geneWindow_ames2nam_local_3allNAM, 
          filename = "ames2nam_3gPCs_allNAM_geneWindow.txt")

# -----------------------------------------
# Import all NAM files
# Convert vcf files to ped files for Gemma
# -----------------------------------------

# Paths of vcf files
chr1.vcf <- "/workdir/mbb262/genotypes/nam/imputed/combined_ibm_nam_all_imputed_bcftools/filtered_DR2_0.8/nam_ibm_imputed_filteredDR2_0.8_chr1.vcf.gz"
chr2.vcf <- "/workdir/mbb262/genotypes/nam/imputed/combined_ibm_nam_all_imputed_bcftools/filtered_DR2_0.8/nam_ibm_imputed_filteredDR2_0.8_chr2.vcf.gz"
chr3.vcf <- "/workdir/mbb262/genotypes/nam/imputed/combined_ibm_nam_all_imputed_bcftools/filtered_DR2_0.8/nam_ibm_imputed_filteredDR2_0.8_chr3.vcf.gz"
chr4.vcf <- "/workdir/mbb262/genotypes/nam/imputed/combined_ibm_nam_all_imputed_bcftools/filtered_DR2_0.8/nam_ibm_imputed_filteredDR2_0.8_chr4.vcf.gz"
chr5.vcf <- "/workdir/mbb262/genotypes/nam/imputed/combined_ibm_nam_all_imputed_bcftools/filtered_DR2_0.8/nam_ibm_imputed_filteredDR2_0.8_chr5.vcf.gz"
chr6.vcf <- "/workdir/mbb262/genotypes/nam/imputed/combined_ibm_nam_all_imputed_bcftools/filtered_DR2_0.8/nam_ibm_imputed_filteredDR2_0.8_chr6.vcf.gz"
chr7.vcf <- "/workdir/mbb262/genotypes/nam/imputed/combined_ibm_nam_all_imputed_bcftools/filtered_DR2_0.8/nam_ibm_imputed_filteredDR2_0.8_chr7.vcf.gz"
chr8.vcf <- "/workdir/mbb262/genotypes/nam/imputed/combined_ibm_nam_all_imputed_bcftools/filtered_DR2_0.8/nam_ibm_imputed_filteredDR2_0.8_chr8.vcf.gz"
chr9.vcf <- "/workdir/mbb262/genotypes/nam/imputed/combined_ibm_nam_all_imputed_bcftools/filtered_DR2_0.8/nam_ibm_imputed_filteredDR2_0.8_chr9.vcf.gz"
chr10.vcf <- "/workdir/mbb262/genotypes/nam/imputed/combined_ibm_nam_all_imputed_bcftools/filtered_DR2_0.8/nam_ibm_imputed_filteredDR2_0.8_chr10.vcf.gz"

# vcf to gds
setwd("/workdir/mbb262/data")
all_chrs_vcf <- snpgdsVCF2GDS(c(chr1.vcf, chr2.vcf, chr3.vcf, chr4.vcf, chr5.vcf, chr6.vcf, chr7.vcf, chr8.vcf, chr9.vcf, chr10.vcf), 
                              "all_nam_chrs.gds", 
                              method="copy.num.of.ref", 
                              snpfirstdim=TRUE , 
                              ignore.chr.prefix = "S")

# GDS to bed for gemma
snpgdsGDS2BED(all_chrs_vcf, bed.fn="nam_ibm_imputed_filteredDR2_0.8_all_chrs")

# Load in fam file 
setwd("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/data/genos/hapmap_v321_snps")
fam <- read.table("nam_ibm_imputed_filteredDR2_0.8_all_chrs.fam", sep = "\t")
colnames(fam) <- c("V1", "IBM_compat", "V3", "V4", "V5", "V6")

# CAREFULLY merge in NAM pehnotypes and rearrange them
library(plyr)
temp <- plyr::join(fam, all_NAM_phenos, type = "left", match = "first")
fam2 <- temp[,c(2,1,3:5,9:17)]

# Save column names
column_names_fam <- colnames(fam2)

# Export temp in gemma format
write.table(fam2, "wallace_phenotypes_fam_2020_02_25.fam", quote = F, sep = "\t",
            col.names = F, row.names = F)

# Load in 3 gPC covariates and format
gpc <- read.delim("ames2nam_3gPCs_allNAM_gPCs_only_R.txt")
colnames(gpc)[1] <- "IBM_compat"
gpc_fam <- plyr::join(fam, gpc, type = "left", match = "first")
gpc_fam <- cbind(rep(1, nrow(gpc_fam)) ,gpc_fam[,c(7:9)])

write.table(gpc_fam, "3gPCs_for_gemma.txt", quote = F, sep = "\t",
            col.names = F, row.names = F)


# ---------------------------------------
# Plot results from the following models
# ---------------------------------------

# Set directories
setwd("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/results/2020_02_24_newHapmapSNPs")
setwd("/workdir/mbb262/results/2020_02_24_newHapmapSNPs")

# 1- y ~ SNP + 3 gPCs
model1_chr1 <- read.delim("model1_chrom1_ames2nam_3gPCs_allNAM_gPCs_FAresults_2020_02_24.txt", header = T)
model1_chr2 <- read.delim("model1_chrom2_ames2nam_3gPCs_allNAM_gPCs_FAresults_2020_02_24.txt", header = T)
model1_chr3 <- read.delim("model1_chrom3_ames2nam_3gPCs_allNAM_gPCs_FAresults_2020_02_24.txt", header = T)
model1_chr4 <- read.delim("model1_chrom4_ames2nam_3gPCs_allNAM_gPCs_FAresults_2020_02_24.txt", header = T)
model1_chr5 <- read.delim("model1_chrom5_ames2nam_3gPCs_allNAM_gPCs_FAresults_2020_02_24.txt", header = T)
model1_chr6 <- read.delim("model1_chrom6_ames2nam_3gPCs_allNAM_gPCs_FAresults_2020_02_24.txt", header = T)
model1_chr7 <- read.delim("model1_chrom7_ames2nam_3gPCs_allNAM_gPCs_FAresults_2020_02_24.txt", header = T)
model1_chr8 <- read.delim("model1_chrom8_ames2nam_3gPCs_allNAM_gPCs_FAresults_2020_02_24.txt", header = T)
model1_chr9 <- read.delim("model1_chrom9_ames2nam_3gPCs_allNAM_gPCs_FAresults_2020_02_24.txt", header = T)
model1_chr10 <- read.delim("model1_chrom10_ames2nam_3gPCs_allNAM_gPCs_FAresults_2020_02_24.txt", header = T)
model1 <- rbind(model1_chr1, model1_chr2, model1_chr3, model1_chr4, model1_chr5, model1_chr6,
                model1_chr7, model1_chr8, model1_chr9, model1_chr10)
model1_filtered <- model1[which(model1$p < 0.001),]


# 2- y ~ SNP + 3 gPCs + kinship matrix (run in Gemma)
trait1_gemma <- read.delim("~/tassel_gwa/results/2020_03_07_gemma/trait1/output/NAM_gemma_3gPCs_kinship_trait1.assoc.txt")
trait2_gemma <- read.delim("~/tassel_gwa/results/2020_03_07_gemma/trait2/output/NAM_gemma_3gPCs_kinship_trait2.assoc.txt")
trait3_gemma <- read.delim("~/tassel_gwa/results/2020_03_07_gemma/trait3/output/NAM_gemma_3gPCs_kinship_trait3.assoc.txt")
trait4_gemma <- read.delim("~/tassel_gwa/results/2020_03_07_gemma/trait4/output/NAM_gemma_3gPCs_kinship_trait4.assoc.txt")
trait5_gemma <- read.delim("~/tassel_gwa/results/2020_03_07_gemma/trait5/output/NAM_gemma_3gPCs_kinship_trait5.assoc.txt")
trait6_gemma <- read.delim("~/tassel_gwa/results/2020_03_07_gemma/trait6/output/NAM_gemma_3gPCs_kinship_trait6.assoc.txt")
trait7_gemma <- read.delim("~/tassel_gwa/results/2020_03_07_gemma/trait7/output/NAM_gemma_3gPCs_kinship_trait7.assoc.txt")
trait8_gemma <- read.delim("~/tassel_gwa/results/2020_03_07_gemma/trait8/output/NAM_gemma_3gPCs_kinship_trait8.assoc.txt")
trait9_gemma <- read.delim("~/tassel_gwa/results/2020_03_07_gemma/trait9/output/NAM_gemma_3gPCs_kinship_trait9.assoc.txt")

trait8_gemma$annotation

# 3- y ~ SNP + 3 gPCs + 660 wPCs
model3_chr1 <- read.delim("model3_chrom1_ames2nam_3gPCs_allNAM_geneWindow_FAresults_2020_02_24.txt", header = T)
model3_chr2 <- read.delim("model3_chrom2_ames2nam_3gPCs_allNAM_geneWindow_FAresults_2020_02_24.txt", header = T)
model3_chr3 <- read.delim("model3_chrom3_ames2nam_3gPCs_allNAM_geneWindow_FAresults_2020_02_24.txt", header = T)
model3_chr4 <- read.delim("model3_chrom4_ames2nam_3gPCs_allNAM_geneWindow_FAresults_2020_02_24.txt", header = T)
model3_chr5 <- read.delim("model3_chrom5_ames2nam_3gPCs_allNAM_geneWindow_FAresults_2020_02_24.txt", header = T)
model3_chr6 <- read.delim("model3_chrom6_ames2nam_3gPCs_allNAM_geneWindow_FAresults_2020_02_24.txt", header = T)
model3_chr7 <- read.delim("model3_chrom7_ames2nam_3gPCs_allNAM_geneWindow_FAresults_2020_02_24.txt", header = T)
model3_chr8 <- read.delim("model3_chrom8_ames2nam_3gPCs_allNAM_geneWindow_FAresults_2020_02_24.txt", header = T)
model3_chr9 <- read.delim("model3_chrom9_ames2nam_3gPCs_allNAM_geneWindow_FAresults_2020_02_24.txt", header = T)
model3_chr10 <- read.delim("model3_chrom10_ames2nam_3gPCs_allNAM_geneWindow_FAresults_2020_02_24.txt", header = T)
model3 <- rbind(model3_chr1, model3_chr2, model3_chr3, model3_chr4, model3_chr5, model3_chr6,
                model3_chr7, model3_chr8, model3_chr9, model3_chr10)
model3_filtered <- model3[which(model3$p < 0.001),]

# Save a test batch of SNPs to find snp effects
# test_model3 <- model3[which(model3$p < 0.00000001),] # -log10(0.00000001) = 8
# write.table(test_model3, "test_model3_results_8_fastAssociation.txt", row.names = F)
library(dplyr)
dts_model3 <- model3 %>% filter(Trait == "Days_To_Silk_BLUP_Sum0607_Buckler2009")
write.csv(dts_model3, "dts_model3.csv")

# 4- y ~ SNP + 26 gPCs + 660 wPCs
model4_chr1 <- read.delim("model4_chrom1_ames2nam_26gPCs_allNAM_geneWindow_FAresults_2020_02_24.txt", header = T)
model4_chr2 <- read.delim("model4_chrom2_ames2nam_26gPCs_allNAM_geneWindow_FAresults_2020_02_24.txt", header = T)
model4_chr3 <- read.delim("model4_chrom3_ames2nam_26gPCs_allNAM_geneWindow_FAresults_2020_02_24.txt", header = T)
model4_chr4 <- read.delim("model4_chrom4_ames2nam_26gPCs_allNAM_geneWindow_FAresults_2020_02_24.txt", header = T)
model4_chr5 <- read.delim("model4_chrom5_ames2nam_26gPCs_allNAM_geneWindow_FAresults_2020_02_24.txt", header = T)
model4_chr6 <- read.delim("model4_chrom6_ames2nam_26gPCs_allNAM_geneWindow_FAresults_2020_02_24.txt", header = T)
model4_chr7 <- read.delim("model4_chrom7_ames2nam_26gPCs_allNAM_geneWindow_FAresults_2020_02_24.txt", header = T)
model4_chr8 <- read.delim("model4_chrom8_ames2nam_26gPCs_allNAM_geneWindow_FAresults_2020_02_24.txt", header = T)
model4_chr9 <- read.delim("model4_chrom9_ames2nam_26gPCs_allNAM_geneWindow_FAresults_2020_02_24.txt", header = T)
model4_chr10 <- read.delim("model4_chrom10_ames2nam_26gPCs_allNAM_geneWindow_FAresults_2020_02_24.txt", header = T)
model4 <- rbind(model4_chr1, model4_chr2, model4_chr3, model4_chr4, model4_chr5, model4_chr6,
                model4_chr7, model4_chr8, model4_chr9, model4_chr10)
model4_filtered <- model4[which(model4$p < 0.001),]


# ---------------------------------------
# Make combined table for loci captured
# Using quantile p-values
# ---------------------------------------

source("~/git_projects/haplotype_v_snp_gwas/src/R/pca_testing/count_gwas_overlaps.R")

# Get prior ranges
# priors <- prior_verified_gene_intervals(gene_ld_buffer = 2500000)
priors_gene <- prior_verified_gene_intervals(gene_ld_buffer = 50000, gene_type = "prior")
priors_random <- prior_verified_gene_intervals(gene_ld_buffer = 50000, gene_type = "random")

# Model 1
model1_overlaps_gene <- count_sig_gwas_overlaps_quantile(gwas_results = model1,
                                                    prior_ranges = priors_gene, max_snp_overlap = 1, df_id = "model1")
model1_overlaps_random <- count_sig_gwas_overlaps_quantile(gwas_results = model1,
                                                    prior_ranges = priors_random, max_snp_overlap = 1, df_id = "model1")

# Model 2

# Model 3
model3_overlaps_gene <- count_sig_gwas_overlaps_quantile(gwas_results = model3, prior_ranges = priors_gene, 
                                                    max_snp_overlap = 1, df_id = "model3")
model3_overlaps_random <- count_sig_gwas_overlaps_quantile(gwas_results = model3, prior_ranges = priors_random, 
                                                    max_snp_overlap = 1, df_id = "model3")

# Model 4
model4_overlaps_gene <- count_sig_gwas_overlaps_quantile(gwas_results = model4, prior_ranges = priors_gene,
                                                    max_snp_overlap = 1, df_id = "model4")
model4_overlaps_random <- count_sig_gwas_overlaps_quantile(gwas_results = model4, prior_ranges = priors_random,
                                                    max_snp_overlap = 1, df_id = "model4")


# Compile and save, DO NOT JUST CBIND ROWS TOGETHER USE MERGE 
#         - omg merrit no wonder why everything looked random before
temp1 <- merge.data.frame(model1_overlaps_gene, model1_overlaps_random, by = "trait", suffixes = c("_priorGene","_randomGene"))
temp3 <- merge.data.frame(model3_overlaps_gene, model3_overlaps_random, by = "trait", suffixes = c("_priorGene","_randomGene"))
temp4 <- merge.data.frame(model4_overlaps_gene, model4_overlaps_random, by = "trait", suffixes = c("_priorGene","_randomGene"))

# TUrn into numeric values
temp1[,2:11] <- apply(temp1[,2:11], 2, function(x) as.numeric(as.character(x)))
temp3[,2:11] <- apply(temp3[,2:11], 2, function(x) as.numeric(as.character(x)))
temp4[,2:11] <- apply(temp4[,2:11], 2, function(x) as.numeric(as.character(x)))

# Get datasetup for chi-squared test
chi_model1 <- temp1[,c(1,2,4,7,9)]
chi_model3 <- temp3[,c(1,2,4,7,9)]
chi_model4 <- temp4[,c(1,2,4,7,9)]

# if value is 0 replace with 1 just for this test
chi_model1[chi_model1==0] <-1
chi_model3[chi_model3==0] <-1
chi_model4[chi_model4==0] <-1

# Calculate chi-squared
df <- chi_model1
holder_m1_1 <- data.frame() # cols 2,4
holder_m1_01 <- data.frame() # cols 3,5

for (trait in seq(1:9)){
  # Do the test
  temp <- chisq.test(df[trait,c(3,5)])$p.value
  
  # Make it a dataframe
  lala <- data.frame(df[trait,1], temp)
  
  # Combine results
  holder_m1_01 <- rbind(holder_m1_01, lala)
}
  

all_chi <- holder_m1_1
all_chi <- cbind(holder_m1_1, holder_m1_01)
all_chi <- all_chi[,c(1,2,4)]
colnames(all_chi) <- c("Trait", "Model1_1", "Model1_01")
write.csv(all_chi, file = "Model3_chi_tests.csv", row.names = F)


# ---------------------------------
#         plotting results
# ---------------------------------

unique_phenos <- unique(model1$Trait)
phenos <- 1
ylim_gwa <- max(-log10(model4_filtered[which(model4_filtered$Trait == unique_phenos[phenos]), ]$p)) + 3

model1_plot <- manhattan_plot(model1_filtered[which(model1_filtered$Trait == unique_phenos[phenos]), ],
                              sig_threshold = 15, ylim = ylim_gwa, title = paste0("Model 1: y ~ SNP + 3gPCs for \n ", unique_phenos[phenos]))
model2_plot <- manhattan_plot(model1_filtered[which(model1_filtered$Trait == unique_phenos[phenos]), ],
                              sig_threshold = 15, ylim = ylim_gwa, title = paste0("Model 2: y ~ SNP + 3gPCs + kinship for \n", unique_phenos[phenos]))
model3_plot <- manhattan_plot(model3_filtered[which(model3_filtered$Trait == unique_phenos[phenos]), ],
                              sig_threshold = 15, ylim = ylim_gwa, title = paste0("Model 3: y ~ SNP + 3gPCs + 660wPCs for \n ", unique_phenos[phenos]))
model4_plot <- manhattan_plot(model4_filtered[which(model4_filtered$Trait == unique_phenos[phenos]), ],
                              sig_threshold = 15, ylim = ylim_gwa, title = paste0("Model 4: y ~ SNP + 26gPCs + 660wPCs for \n", unique_phenos[phenos]))

# Combine plots and save
# (model1_plot + model2_plot)/(model3_plot + model4_plot)
(model1_plot + model3_plot + model4_plot)


# ----------------------
# Plot only DTA and DTS
# ----------------------

# Plot singles
model <- model3

# Subset by trait
model_fum <- model[which(model$Trait == "Fumarate_Blup_Wallace2014" & model$p < 0.5),]
model_dts <- model[which(model$Trait == "Days_To_Silk_BLUP_Sum0607_Buckler2009" & model$p < 0.01),]
model_cob <- model[which(model$Trait == "cob_diameter_raw_Hung2012" & model$p < 0.01),]

# model3_dts <- model3[which(model3$Trait == "Days_To_Silk_BLUP_Sum0607_Buckler2009" & model3$p < 0.000001),]
# model4_dts <- model4[which(model4$Trait == "Days_To_Silk_BLUP_Sum0607_Buckler2009" & model4$p < 0.000001),]

# Manhattan Plot
a_fum <- manhattan_plot(model_fum, sig_threshold_1 = 2.027, sig_threshold_2 = 4.13, ylim = 20, title = "Fumarate")
a_dts <- manhattan_plot(model_dts, sig_threshold_1 = 6.45, sig_threshold_2 = 19.46, ylim = 45, title = "Days to Silk")
a_cob <- manhattan_plot(model_cob, sig_threshold_1 = 3.85, sig_threshold_2 = 8.26, ylim =20, title = "Cob Diameter")

# m3_dts <- manhattan_plot(model3_dts, sig_threshold = 19.46, ylim = 115, title = "DTS Model 3")
# m4_dts <- manhattan_plot(model4_dts, sig_threshold = 20.35, ylim = 115, title = "DTS Model 4")

# Use patchwork
# a_fum+a_dts
# a_dts + m3_dts + m4_dts
# temp <- a_fum/a_dts/a_cob
# ggsave(filename = "complexity_model3_pcs.png", plot = temp, width = 9, units = "in")
ggsave(filename = "fumarate_model3_pcs.png", plot = a_fum, width = 6, units = "in")
ggsave(filename = "DTS_model3_pcs.png", plot = a_dts, width =6, units = "in")
ggsave(filename = "cob_model3_pcs.png", plot = a_cob, width = 6, units = "in")


