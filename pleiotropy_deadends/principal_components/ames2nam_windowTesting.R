# ---------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2020-03-16
#
# Description 
#   - Parameter tuning for GWAS models
#   - Trying different iterations of windows
# 
# ---------------------------------------------------------------

# Load in source scripts
source('~/git_projects/pleiotropy/src/R/03_ames2any_matrix.R')
source('~/git_projects/pleiotropy/src/R/04_local_window_funs.R')
source("~/git_projects/pleiotropy/src/R/pca_testing/count_gwas_overlaps.R")
source('~/git_projects/pleiotropy/src/R/random_scripts/manhattan_plot.R')
source('~/git_projects/pleiotropy/src/R/05_create_gene_windows.R')

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

write.csv(global_ames_BV_3gPCs$B, "global_ames_BV_3gPCs_B_matrix.csv")
write.csv(global_ames_BV_3gPCs$V, "global_ames_BV_3gPCs_V_matrix.csv")

# Transfer loadings and coefficients over to all X (Ames, NAM) to get adjusted global PCs
# Get Q matrix (population structure)
ames2nam_gPCs_3gPCs_allNAM <- ames2any(X_all = X_allNAM$X_all,
                                       Q = NULL, B = global_ames_BV_3gPCs$B, V = global_ames_BV_3gPCs$V)
Q_all_3gPCs_allNAM <- cbind(1, ames2nam_gPCs_3gPCs_allNAM)

write.csv(Q_all_3gPCs_allNAM, "Q_all_3gPCs_allNAM.csv")

# ------------------------------------------
#   Calculate local PCs in defined windows
#           formatting step
# ------------------------------------------

# Use function
window_file <- create_windows(362)

# 3 gPCs, 220 windows = 181 genes/window
window_file181 <- create_windows(181)

# -------------------------------------------
# Local PC function with gene number windows
# -------------------------------------------

# Use on variable gene number windows with ames to nam PCs
# 3 gPCs, 114 windows = 362 genes/window
geneWindow_ames2nam_local_3allNAM <- calc_local_pcs_gene_window(X_all = X_allNAM,
                                                                genofileID_ames = genofile_ames,
                                                                Q_all = Q_all_3gPCs_allNAM,
                                                                numLocalPcs = 3,
                                                                windowFile = window_file)

geneWindow181_ames2nam_local_3allNAM <- calc_local_pcs_gene_window(X_all = X_allNAM,
                                                                   genofileID_ames = genofile_ames,
                                                                   Q_all = Q_all_3gPCs_allNAM,
                                                                   numLocalPcs = 3,
                                                                   windowFile = window_file181)

# 3 gPCs, 220 windows = 181 genes/window, 2 PCs per window

# 3 gPCs, 220 windows = 181 genes/window, 1 PCs per window


# Get phenotypes, tasselize, save
all_NAM_phenos <- read.table("~/git_projects/haplotype_v_snp_gwas/data/all_NAM_phenos.txt", sep="", header = TRUE)
all_NAM_phenos <- all_NAM_phenos[,c(2,4,17,9)]

# "Fix" IBM names
all_NAM_phenos$IBM_compat <- as.character(all_NAM_phenos$Geno_Code)
all_NAM_phenos[c(3203:3402),5] <- as.character(all_NAM_phenos[c(3203:3402),2])
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
  # together <- together[,c(13,4:12,14:ncol(together))]
  print(colnames(together))
  together <- together[,c(6,4:5,7:ncol(together))]
  
  # Change name of first column from taxa ID to trait
  # colnames(together)[1] <- "<Trait>"
  
  # Remove missing and duplicated values
  together <- na.omit(together)
  together <- together[!duplicated(together),]
  
  # To export dataset as is to work on in R
  # return(together)
  
  # Export global PCs only  
  # Subset out only the global PCs
  global <- together[,c(1,4:6)]
  global[,1] <- as.character(global[,1])
  colnames(global)[1] <- "<Trait>"
  header_global <- c("<Covariate>", rep("", ncol(global)-1))
  global <- rbind(header_global, colnames(global), global)
  write.table(global, file = "ames2nam_3gPCs_allNAM_gPCs_only.txt",
              row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  
  # Export phenotypes
  # Subset out phenotypes and format
  phenos <- together[,c(1:3)]
  phenos[,1] <- as.character(phenos[,1])
  colnames(phenos)[1] <- "<Trait>"
  # header_phenos <- c("<Trait>", rep("", ncol(phenos)-1))
  # phenos <- rbind(header_phenos, colnames(phenos), phenos)
  write.table(phenos, file = "wallace_phenos_2020_03_16.txt", row.names = FALSE,
              col.names = T, sep = "\t", quote = FALSE)
  
  # Only export window PCs
  together <- together[,c(1,7:ncol(together))]
  together[,1] <- as.character(together[,1])
  colnames(together)[1] <- "<Trait>"
  header <- c("<Covariate>", rep("", ncol(together)-1))
  together <- rbind(header, colnames(together), together)
  write.table(together, file = filename, row.names = FALSE, 
              col.names = FALSE, sep = "\t", quote = FALSE)
}


# Gene window options
tasselize(X_allNAM$X_all, ames2nam_gPCs_3gPCs_allNAM, geneWindow_ames2nam_local_3allNAM, 
          filename = "ames2nam_3gPCs_allNAM_geneWindow_360geneWindow.txt")

tasselize(X_allNAM$X_all, ames2nam_gPCs_3gPCs_allNAM, geneWindow181_ames2nam_local_3allNAM, 
          filename = "ames2nam_3gPCs_allNAM_geneWindow_181geneWindow.txt")


# -------------------------
#     Visualize results
# -------------------------

# Look at windows
vis_windows_genes(window_file_df = window_file)
ggsave("360_genes_per_window.png")
vis_windows_genes(window_file_df = window_file181)
ggsave("181_genes_per_window.png")

# Read in files
setwd("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/results/2020_03_16_variedWindowSizeChangedSNPs")
model5_360_8 <- read.table("model5_chrom8_ames2nam_3gPCs_360geneWindow_FAresults_2020_03_16.txt", header = T)
model5_360_10 <- read.table("model5_chrom10_ames2nam_3gPCs_360geneWindow_FAresults_2020_03_16.txt", header = T)
model6_181_8 <- read.table("model6_chrom8_ames2nam_3gPCs_181geneWindow_FAresults_2020_03_16.txt", header = T)
model6_181_10 <- read.table("model6_chrom10_ames2nam_3gPCs_181geneWindow_FAresults_2020_03_16.txt", header = T)

# Plot the results within calculated PC windows
source('~/git_projects/haplotype_v_snp_gwas/src/R/pca_testing/create_gene_windows.R')

# Chromosome 8 results
a <- manhattan_with_windows(results = model5_360_8,
                            upper_bound = 40,
                            chr = 8,
                            window_file_df = window_file,
                            title ="y ~ SNP + 3gPCs + 114 windows * 3 PCs per window \n for DTS (362 genes/window)")
b <- manhattan_with_windows(results = model6_181_8,
                            upper_bound = 20,
                            chr = 8,
                            window_file_df = window_file181,
                            title ="y ~ SNP + 3gPCs + 220 windows * 3 PCs per window \n for DTS (181 genes/window))")

# Chromosome 10 results
c <- manhattan_with_windows10(results = model5_360_10, 
                              upper_bound = 75, 
                              chr = 10, 
                              window_file_df = window_file, 
                              title ="y ~ SNP + 3gPCs + 114 windows * 3 PCs per window \n for DTS (362 genes/window)")

d <- manhattan_with_windows10(results = model6_181_10, 
                              upper_bound = 75, 
                              chr = 10, 
                              window_file_df = window_file181, 
                              title ="y ~ SNP + 3gPCs + 220 windows * 3 PCs per window \n for DTS (181 genes/window))")

library(patchwork)

(b / a) | (d / c)
ggsave("two models two chromosomes.png")


# ----------------------------------------
# Capture variable number of PCs within
# a window
# ----------------------------------------

# Calculate 3 global PCs for Ames only

# Calculate local PCs for ames only after adjusting by global PCs
# Function to get Ames B and V matrices given original X matrix
get_ames_BV_and_var <- function(X_ames, Q = NULL, num_PCs){
  
  # If X_ames is a "SNPGDSFile" object, generate an X matrix
  if (class(X_ames)[1] == "SNPGDSFileClass"){
    # Turn genotypes into numeric values
    X_list <- snpgdsGetGeno(X_ames, snpfirstdim = FALSE, with.id = TRUE, verbose = FALSE)
    
    # Format X matrix get alternate counts instead reference counts
    X <- X_list$genotype
    rownames(X) <- X_list$sample.id
    colnames(X) <- X_list$snp.id
    X <- 2-X
    
  } else {
    # If X_ames is a matrix already, keep it and rename it
    X <- X_ames
  }
  
  # If Q matrix (population structure) is not present, create matrix of 1's
  if (is.null(Q)) {
    Q <- matrix(1, nrow=nrow(X), ncol=1)
  }
  
  # Get coefficients from effect of Q (pop. str) on X (markers)
  B <- lm(X ~ Q - 1)$coefficients
  
  # Adjust for covariates
  X <- X - Q %*% B
  
  # Do PCA, return rotations
  pca_results <- prcomp(X, center = FALSE, rank. = num_PCs)
  
  V <- pca_results$rotation
  
  # Return B and V matrices
  return(list(B = B, V = V, pca = data.frame(summary(pca_results)[6])))
}

global_ames_BV_var_3gPCs <- get_ames_BV_and_var(X_ames = genofile_ames, Q = NULL, num_PCs = 3)

# Format X matrix get alternate counts instead reference counts
Xames_getgeno <- snpgdsGetGeno(genofile_ames, snpfirstdim = FALSE, with.id = TRUE, verbose = FALSE)
Xames <- Xames_getgeno$genotype
rownames(Xames) <- Xames_getgeno$sample.id
colnames(Xames) <- Xames_getgeno$snp.id
Xames <- 2-Xames
info <- snpgdsSNPList(genofile_ames)
Xames <- list(X_all = Xames, snp_info = info)


# Calculate ames PCS in windows
ames2ames_gPCs_3gPCs <- ames2any(X_all = Xames$X_all, Q = NULL, 
                                 B = global_ames_BV_var_3gPCs$B, V = global_ames_BV_var_3gPCs$V)
Q_all_3gPCs_ames <- cbind(1,ames2ames_gPCs_3gPCs)

write.csv(global_ames_BV_var_3gPCs$B, "global_ames_BV_var_3gPCs_B_matrix_amesOnly.csv")
write.csv(global_ames_BV_var_3gPCs$V, "global_ames_BV_var_3gPCs_V_matrix_amesOnly.csv")
write.csv(Q_all_3gPCs_ames, "Q_all_3gPCs_ames.csv")

# Use the function with new and improved window coordinates
window_file181 <- getWindowCoords(genesPerWindow = 181)
geneWindow181_ames_window <- calc_local_pcs_gene_window_var(X_all = Xames,
                                                        genofileID_ames = genofile_ames,
                                                        Q_all =Q_all_3gPCs_ames,
                                                        numLocalPcs = 10,
                                                        windowFile = window_file181)

window_file360 <- getWindowCoords(genesPerWindow = 360)
geneWindow360_ames_window <- calc_local_pcs_gene_window_var(X_all = Xames,
                                                            genofileID_ames = genofile_ames,
                                                            Q_all =Q_all_3gPCs_ames,
                                                            numLocalPcs = 10,
                                                            windowFile = window_file360)


# -------------------
#      Plotting
# --------------------
# Stacked barchart of variance explained by each Ames PC
library(ggplot2)
a <- ggplot(geneWindow360_ames_window, aes(x = window, y = Proportion.of.Variance, fill = PCID)) + 
  geom_bar(stat = 'identity') +
  facet_grid(cols = vars(chr), scales = "free_x")

# See how many snps are in each window
snps_per_window <- geneWindow360_ames_window %>% filter(PCID == "importance.PC4")
b <- ggplot(snps_per_window, aes(x = window, y = snps_per_window)) + 
  geom_bar(stat = 'identity') +
  facet_grid(cols = vars(chr), scales = "free_x")
a/b

# Scree plot
ggplot(geneWindow181_ames_window, aes(x = PCID, y = Proportion.of.Variance)) +
  geom_boxplot()


# ---------------------
# Dynamic PC selection
# ---------------------

# Select enough windows to get to 15% varaince explained per window

# 181 gene method
temp_window <- lapply(X = seq_len(nrow(window_file181)), FUN = function(i) {
  all_pcs_1_window <- geneWindow181_ames_window %>% filter(windowID == window_file181[i,1])
  # Select enough PCs to get to 15% varaince explained
  filtered <- all_pcs_1_window %>% filter(Cumulative.Proportion >= 0.15)
  filtered <- data.frame(filtered[which.min(filtered$Cumulative.Proportion),]) # Minimum number of rows to get to 15%
  filtered <- filtered[,c(4:5)]
  filtered$PCID <- gsub("importance.PC", "", filtered$PCID)
  filtered$num_pcs <- as.numeric(filtered$PCID)
  return(filtered)
}) %>% do.call("rbind", .)

window_file181 <- merge(x = window_file181, y = temp_window)
window_file181 <- window_file181 %>% arrange(chrom, start)

# 360 gene method
temp_window360 <- lapply(X = seq_len(nrow(window_file360)), FUN = function(i) {
  all_pcs_1_window <- geneWindow360_ames_window %>% filter(windowID == window_file360[i,1])
  # Select enough PCs to get to 15% varaince explained
  filtered <- all_pcs_1_window %>% filter(Cumulative.Proportion >= 0.15)
  filtered <- data.frame(filtered[which.min(filtered$Cumulative.Proportion),]) # Minimum number of rows to get to 15%
  filtered <- filtered[,c(4:5)]
  filtered$PCID <- gsub("importance.PC", "", filtered$PCID)
  filtered$num_pcs <- as.numeric(filtered$PCID)
  return(filtered)
}) %>% do.call("rbind", .)

window_file360 <- merge(x = window_file360, y = temp_window360)
window_file360 <- window_file360 %>% arrange(chrom, start)


# Now calculate variable number of local PCs within NAM using information
# on how many PCs to use per window in Ames (from above)
source('~/git_projects/haplotype_v_snp_gwas/src/R/pca_testing/local_window_funs.R')

geneWindow181_ames2nam_local_3allNAM_varaible <- calc_local_pcs_gene_window_variableNumPCs(X_all = X_allNAM,
                                                                genofileID_ames = genofile_ames,
                                                                Q_all = Q_all_3gPCs_allNAM,
                                                                windowFile = window_file181)

geneWindow360_ames2nam_local_3allNAM_varaible <- calc_local_pcs_gene_window_variableNumPCs(X_all = X_allNAM,
                                                                                        genofileID_ames = genofile_ames,
                                                                                        Q_all = Q_all_3gPCs_allNAM,
                                                                                        windowFile = window_file360)

# Format and Export
temp <- geneWindow360_ames2nam_local_3allNAM_varaible
temp <- data.frame(temp)
temp$taxa <- as.character(rownames(temp))
temp <- temp[,c(ncol(temp), 1:(ncol(temp)-1))]
colnames(temp)[1] <- "<Trait>"
header <- c("<Covariate>", rep("", ncol(temp)-1))
temp <- rbind(header, colnames(temp), temp)
write.table(temp, file = "geneWindow360_ames2nam_local_3allNAM_varaibleNumPCs.txt", 
            row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

# Load in tassel results
m7_chr8 <- read.table("model7_chrom8_ames2nam_3gPCs_181geneWindow_FAresults_2020_03_24.txt", header = TRUE)
m7_chr10 <- read.table("model7_chrom10_ames2nam_3gPCs_181geneWindow_FAresults_2020_03_24.txt", header = TRUE)
m8_chr8 <- read.table("model8_chrom8_ames2nam_3gPCs_360geneWindow_FAresults_2020_03_24.txt", header = TRUE)
m8_chr10 <- read.table("model8_chrom10_ames2nam_3gPCs_360geneWindow_FAresults_2020_03_24.txt", header = TRUE)

# Visualize

# Chromosome 8 results
a <- manhattan_with_windows8(results = m7_chr8,
                            upper_bound = 25,
                            chr = 8,
                            window_file_df = window_file181,
                            title ="y ~ SNP + 3gPCs + 220 windows * 1-6 PCs per window \n for DTS (181 genes/window)")
b <- manhattan_with_windows10(results = m7_chr10,
                            upper_bound = 70,
                            chr = 10,
                            window_file_df = window_file181,
                            title ="y ~ SNP + 3gPCs + 220 windows * 1-6 PCs per window \n for DTS (181 genes/window))")

c <- manhattan_with_windows8(results = m8_chr8,
                        upper_bound = 25,
                        chr = 8,
                        window_file_df = window_file360,
                        title ="y ~ SNP + 3gPCs + 114 windows * 1-8 PCs per window \n for DTS (360 genes/window)")
d <- manhattan_with_windows10(results = m8_chr10,
                         upper_bound = 70,
                         chr = 10,
                         window_file_df = window_file360,
                         title ="y ~ SNP + 3gPCs + 114 windows * 1-8 PCs per window \n for DTS (360 genes/window))")

(c/a) | (d/b)
ggsave("equal varaince in wPCs and modified window function.png")


# ----------------------------------------------
# Calculate PCs with additional midway windows
# ----------------------------------------------

# Create window file with midway windows
source('~/git_projects/haplotype_v_snp_gwas/src/R/pca_testing/create_gene_windows.R')
window_file_midway_360 <- getWindowCoords_and_midwayWindows(genesPerWindow = 360)

# Calculate window PCs within Ames, see how many I need for NAM
source('~/git_projects/haplotype_v_snp_gwas/src/R/pca_testing/local_window_funs.R')
geneWindow_midway_360_ames_window <- calc_local_pcs_gene_window_var(X_all = Xames,
                                                            genofileID_ames = genofile_ames,
                                                            Q_all =Q_all_3gPCs_ames,
                                                            numLocalPcs = 10,
                                                            windowFile = window_file_midway_360)

# Plot ames PC plots
# a <- ggplot(geneWindow_midway_360_ames_window, aes(x = window, y = Proportion.of.Variance, fill = PCID)) + 
#   geom_bar(stat = 'identity') +
#   facet_grid(cols = vars(chr), scales = "free_x")
# 
# # Scree plot
# b <- ggplot(geneWindow_midway_360_ames_window, aes(x = PCID, y = Proportion.of.Variance)) +
#   geom_boxplot()
# a/b
# ggsave("scree and variance midway windows.png")

# Pick enough PCs to explain 15% of the variance
midway_window360 <- lapply(X = seq_len(nrow(window_file_midway_360)), FUN = function(i) {
  all_pcs_1_window <- geneWindow_midway_360_ames_window %>% filter(windowID == window_file_midway_360[i,1])
  # Select enough PCs to get to 15% varaince explained
  filtered <- all_pcs_1_window %>% filter(Cumulative.Proportion >= 0.15)
  filtered <- data.frame(filtered[which.min(filtered$Cumulative.Proportion),]) # Minimum number of rows to get to 15%
  filtered <- filtered[,c(4:5)]
  filtered$PCID <- gsub("importance.PC", "", filtered$PCID)
  filtered$num_pcs <- as.numeric(filtered$PCID)
  return(filtered)
}) %>% do.call("rbind", .)

# Merge, arrange, and export for future analyses
window_file_midway_360 <- merge(x = window_file_midway_360, y = midway_window360)
window_file_midway_360 <- window_file_midway_360 %>% arrange(chrom, start)
write.csv(window_file_midway_360, "window_file_midway_360_2.csv", row.names = FALSE)

# Subset
main_window_file_midway_360 <- window_file_midway_360 %>% filter(windowLocation == "mainWindow")
midway_window_file_midway_360 <- window_file_midway_360 %>% filter(windowLocation == "midwayWindow")

# Calculte dynamic number of PCs within NAM
main_geneWindow360_ames2nam_local_3allNAM_varaible <- calc_local_pcs_gene_window_variableNumPCs(X_all = X_allNAM,
                                                                                           genofileID_ames = genofile_ames,
                                                                                           Q_all = Q_all_3gPCs_allNAM,
                                                                                           windowFile = main_window_file_midway_360)
midway_geneWindow360_ames2nam_local_3allNAM_varaible <- calc_local_pcs_gene_window_variableNumPCs(X_all = X_allNAM,
                                                                                           genofileID_ames = genofile_ames,
                                                                                           Q_all = Q_all_3gPCs_allNAM,
                                                                                           windowFile = midway_window_file_midway_360)
# Check window IDs between files, midway should have higher numbers
str(main_geneWindow360_ames2nam_local_3allNAM_varaible)
str(midway_geneWindow360_ames2nam_local_3allNAM_varaible)

# Save in csv format
temp <- midway_geneWindow360_ames2nam_local_3allNAM_varaible
temp <- data.frame(temp)
write.csv(temp, "midwayGeneWindow360_ames2nam_local_3allNAM_varaibleNumPCs.csv")


# Save and export to tassel
temp <- main_geneWindow360_ames2nam_local_3allNAM_varaible
temp <- data.frame(temp)
temp$taxa <- as.character(rownames(temp))
temp <- temp[,c(ncol(temp), 1:(ncol(temp)-1))]
colnames(temp)[1] <- "<Trait>"
header <- c("<Covariate>", rep("", ncol(temp)-1))
temp <- rbind(header, colnames(temp), temp)
write.table(temp, file = "mainGeneWindow360_ames2nam_local_3allNAM_varaibleNumPCs.txt", 
            row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

# Import results
setwd("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/results/2020_03_28_mainAndMidwayWindows")
main_m10_chr8 <- read.table("model10_chrom8_ames2nam_3gPCs_main360geneWindow_FAresults_2020_03_28.txt", header = TRUE)
main_m10_chr10 <- read.table("model10_chrom10_ames2nam_3gPCs_main360geneWindow_FAresults_2020_03_28.txt", header = TRUE)

midway_m10_chr8 <- read.table("model10_chrom8_ames2nam_3gPCs_midway360geneWindow_FAresults_2020_03_28.txt", header = TRUE)
midway_m10_chr10 <- read.table("model10_chrom10_ames2nam_3gPCs_midway360geneWindow_FAresults_2020_03_28.txt", header = TRUE)

# Subset out only DTS
main_m10_chr8 <- main_m10_chr8 %>% filter(Trait == "Days_To_Silk_BLUP_Sum0607_Buckler2009")
main_m10_chr10 <- main_m10_chr10 %>% filter(Trait == "Days_To_Silk_BLUP_Sum0607_Buckler2009")
midway_m10_chr8 <- midway_m10_chr8 %>% filter(Trait == "Days_To_Silk_BLUP_Sum0607_Buckler2009")
midway_m10_chr10 <- midway_m10_chr10 %>% filter(Trait == "Days_To_Silk_BLUP_Sum0607_Buckler2009")

# Merge datasets by factorized position
main_midway_chr8 <- merge(main_m10_chr8, midway_m10_chr8, by = "Marker", suffixes = c("_main", "_midway"))
main_midway_chr10 <- merge(main_m10_chr10, midway_m10_chr10, by = "Marker", suffixes = c("_main", "_midway"))

# Find min p-value row (out of two columns) and save a new column with lowest p-value
main_midway_chr8$p <- apply(main_midway_chr8[,c(7,13)],1,min,na.rm=TRUE)
main_midway_chr10$p <- apply(main_midway_chr10[,c(7,13)],1,min,na.rm=TRUE)

# Format windowfiles to plot differently

# Plot
source('~/git_projects/haplotype_v_snp_gwas/src/R/pca_testing/create_gene_windows.R')
e <- manhattan_with_windows8(results = main_midway_chr8,
                        upper_bound = 25,
                        chr = 8,
                        window_file_df = window_file_midway_360,
                        title ="y ~ SNP + 3gPCs + 218 overlapping windows * 1-8 PCs per window \n for DTS")
f <- manhattan_with_windows10(results = main_midway_chr10,
                         upper_bound = 75,
                         chr = 10,
                         window_file_df = window_file_midway_360,
                         title ="y ~ SNP + 3gPCs + 218 overlapping windows * 1-8 PCs per window \n for DTS")

# (c/e) | (d/f)
e/f
ggsave("equal varaince in wPCs and overlapping wPCs in 10cm models.png")

