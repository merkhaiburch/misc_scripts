# ---------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2020-06-03 
#
# Description 
#   - Gather vcf files with 282 Hapmap3.2.1 genotypes
# - filter to common set of SNPs between Ames, NAM, and others
# - Merge files, export, use to calculate PCs in R
# ---------------------------------------------------------------


# Install package
if (!require("devtools")) install.packages("devtools")
devtools::install_bitbucket(
    repo = "bucklerlab/rTASSEL",
    ref = "master",
    build_vignettes = FALSE
)

# Setting memory
options(java.parameters = c("-Xmx220g"))
temp <- options() # Check to see if memory was set
temp$java.parameters

# Load package
library(rTASSEL)

# Start logging
rTASSEL::startLogger(fullPath = "/workdir/mbb262/goodman282/filtered_for_pca", fileName = NULL)

# Call needed function
rJC <- rJava::J("net/maizegenetics/plugindef/GenerateRCode")

# Set working directory
setwd("/workdir/mbb262/goodman282/filtered_for_pca")

# Load in source scripts
source('~/git_projects/haplotype_v_snp_gwas/src/R/pca_testing/ames2any_matrix.R')
source('~/git_projects/haplotype_v_snp_gwas/src/R/random_scripts/manhattan_plot.R')
source("~/git_projects/haplotype_v_snp_gwas/src/R/pca_testing/create_gene_windows.R")

# Load packages
library(gdsfmt)
library(SNPRelate)
library(dplyr)
library(ggplot2)

# --------------
# Get data
# --------------

# vcf files from blfs1 
#   (/data1/users/mbb262/genotypes/goodman282)
# Ames sites in a position list file 
#   (/data1/users/mbb262/genotypes/ames/imputed/ames_for_pca/ames_for_pca_merged_position_list.json.gz)


# Function to map traits
map_trait <- function(pc_df, num_covars, id_paste = "") {
  
  # Load in phenotypes
  ft <- read.table("/workdir/mbb262/goodman282/flowering_buckler_282pop.txt", sep = "\t", header = TRUE)
  
  # Merge phenotypes with covariates
  ft_pcs_ames2goodman <- merge(x = ft, y = pc_df, by.x = "hapmap_names", by.y = "Taxa")
  colnames(ft_pcs_ames2goodman)[1] <- "Taxa"
  
  # Remove subpopulation
  ft_pcs_ames2goodman_2 <- ft_pcs_ames2goodman %>% select(-Subpopulation, -)
  print(head(ft_pcs_ames2goodman_2))
  
  # Tasselize merged phenotypes + global PCs + window PCs
  tasPhenoDF <- rTASSEL::readPhenotypeFromDataFrame(
    phenotypeDF = ft_pcs_ames2goodman_2,
    taxaID = "Taxa",
    attributeTypes = c(rep("data", 1), rep("covariate", num_covars)))
  message("I am here")
  # For loop to iterate through genotypes/chromosomes
  for (i in seq_len(10)){
    message("I am on chromosome ", i)
    
    # Load in genotype table
    vcf <-  rTASSEL::readGenotypeTableFromPath(path = paste("/workdir/mbb262/rtassel/goodman282/","hmp321_282_agpv4_merged_chr", i, "_imputed_goodman282.vcf.gz", sep = ""),
                                               keepDepth = FALSE)
    print(vcf)
    
    # Join genotypes with (phenotypes + g PCs)
    tasPhenoDF <- rTASSEL::readGenotypePhenotype(
      genoPathOrObj = vcf,
      phenoPathDFOrObj = tasPhenoDF)
    print(tasPhenoDF)
    
    # Do a light MAF filter to remove invariant sites
    tasGenoPhenoFilt <- rTASSEL::filterGenotypeTableSites(
      tasObj = tasPhenoDF,
      siteRangeFilterType = "none",
      siteMinAlleleFreq = 0.01,
      siteMaxAlleleFreq = 1.0)
    
    # Run fast association, write files to disk
    rJC$fastAssociation(
      tasGenoPhenoFilt@jTasselObj,             # TASSEL object
      rJava::.jnew("java.lang.Double", "0.01"), # maxP
      rJava::.jnew("java.lang.Integer", 40L),   # maxThreads
      TRUE,                                    # writeToFile
      paste("chrom_", i, id_paste, sep = "")) # outputFile
  }
  
}


# ------------------------------------------
#     Calculate ames2goodman global PCs
#     uses ames2any_matrix functions
# ------------------------------------------

# Load in Beagle imputed SNPs that were filtered with bcftools
ames.vcf <- "ames_for_pca_merged.vcf"
ames_vcf <- snpgdsVCF2GDS(ames.vcf, "ames.gds", method="copy.num.of.ref", snpfirstdim=TRUE, ignore.chr.prefix = "S")
snpgdsSummary("ames.gds")
genofile_ames <- snpgdsOpen(ames_vcf)

# Load in 282 SNPs that are shared (intersect) with the Ames SNPs
goodman.vcf <- "goodman282_by_ames_sites_for_pca_allChroms.vcf"
goodman_vcf <- snpgdsVCF2GDS(goodman.vcf, "goodman.gds", method="copy.num.of.ref", snpfirstdim=TRUE, ignore.chr.prefix = "S")
snpgdsSummary("goodman.gds")
genofile_goodman <- snpgdsOpen(goodman_vcf)

# Combine genofile objects
X_all <- combine_genofile_matrix(genofile_ames = genofile_ames, genofile_any = genofile_goodman)

# Calculate Ames global PCs, get loadings
global_ames_BV_3gPCs <- get_ames_BV(X_ames = genofile_ames, Q = NULL, num_PCs = 5)

# Transfer loadings and coefficients over to all X (Ames, goodman) to get adjusted global PCs
# Get Q matrix (population structure)
ames2goodman_gPCs_3gPCs_allNAM <- ames2any(X_all = X_all$X_all,
                                           Q = NULL, B = global_ames_BV_3gPCs$B, V = global_ames_BV_3gPCs$V)
write.csv(ames2goodman_gPCs_3gPCs_allNAM, "Q_all_3gPCs_allGoodman.csv")
write.table(ames2goodman_gPCs_3gPCs_allNAM, "Q_all_3gPCs_allGoodman.txt", quote = F, row.names = T, sep = "\t")


# Check varaince
# lala <- get_ames_BV_all_pcs(genofile_ames, num_PCs = 5)
# summary(lala$pca)

# Load in pre-analyzed PCs
ames2goodman_gPCs_3gPCs_allNAM_df <- read.csv("Q_all_3gPCs_allGoodman.csv", header = TRUE) %>% select(-PC4, -PC5)
colnames(ames2goodman_gPCs_3gPCs_allNAM_df) <- c("Taxa", "PC1", "PC2", "PC3")

# ------------------
# Visualize PCs
# ------------------

# Load in data
pcs <- read.csv("~/Downloads/Q_all_3gPCs_allGoodman.csv", header = TRUE)
ft <- read.table("~/Downloads/flowering_buckler_282pop.txt", sep = "\t", header = TRUE)

# Combine data
pcs_df <- merge(x = pcs, y = ft, by.x = "X", by.y = "X.Trait.")

# Packages
library(ggplot2)
library(patchwork)

# Plot
a <- ggplot(pcs_df, aes(x = PC1, y = PC2, color = Subpopulation)) +
  geom_point() +
  ggtitle("PC1 vs PC2: SS vs non-SS") +
  theme(legend.position="bottom")

b <- ggplot(pcs_df, aes(x = PC2, y = PC3, color = Subpopulation)) +
  geom_point() +
  ggtitle("PC2 vs PC3: tropical vs sweet") +
  theme(legend.position="bottom")

c <- ggplot(pcs_df, aes(x = PC1, y = PC3, color = Subpopulation)) +
  geom_point() +
  ggtitle("PC1 vs PC3: ss vs everything") +
  theme(legend.position="bottom")

a|b|c
ggsave("ames2goodman_pcs.png")

# ------------------------------
# Ames to goodman method, 3 PCs
# ------------------------------

map_trait(pc_df = ames2goodman_gPCs_3gPCs_allNAM_df, num_covars = 3, 
          id_paste = "_fast_assoc_results_ames2goodman_flowering_time_buckler_282")

# Load in results
my_files <- list.files(pattern = "fast_assoc_results_ames2goodman_flowering_time_buckler_282.txt$")
dts_ames2goodman_all_dts <- lapply(my_files, data.table::fread) %>% data.table::rbindlist()
dts_ames2goodman_3pcs <- data.frame(dts_ames2goodman_all_dts)

# Get significance thresholds
sig_threshold_1 = -log10(quantile(dts_ames2goodman_3pcs$p, probs = 0.05))
sig_threshold_2 = -log10(quantile(dts_ames2goodman_3pcs$p, probs = 0.01))
ames2goodman_gwas_3pcs <- manhattan_plot(gwas_result = dts_ames2goodman_3pcs, 
                                         sig_threshold_1 = sig_threshold_1 , sig_threshold_2 = sig_threshold_2,
                                         ylim = 25,title = "DTS ~ SNP + 3 PCs (ames to 282)")
# ggsave("DTS_ames2goodman.png")


# ------------------------------
# Ames to goodman method, 2 PCs
# ------------------------------

# Subset down to only 2 PCs
ames2goodman_gPCs_2gPCs_allNAM_df <- read.csv("Q_all_3gPCs_allGoodman.csv", header = TRUE) %>% select(-PC3, -PC4, -PC5)
head(ames2goodman_gPCs_2gPCs_allNAM_df)
colnames(ames2goodman_gPCs_2gPCs_allNAM_df) <- c("Taxa", "PC1", "PC2")


map_trait(pc_df = ames2goodman_gPCs_2gPCs_allNAM_df, num_covars = 2, 
          id_paste = "_fast_assoc_results_ames2goodman_flowering_time_buckler_282_2pcs")

# Load in results
my_files <- list.files(pattern = "_fast_assoc_results_ames2goodman_flowering_time_buckler_282_2pcs.txt$")
my_files
dts_ames2goodman_all_dts <- lapply(my_files, data.table::fread) %>% data.table::rbindlist()
dts_ames2goodman_2pcs <- data.frame(dts_ames2goodman_all_dts)

# Get significance thresholds
sig_threshold_1 = -log10(quantile(dts_ames2goodman_2pcs$p, probs = 0.05))
sig_threshold_2 = -log10(quantile(dts_ames2goodman_2pcs$p, probs = 0.01))
ames2goodman_gwas_2pcs <- manhattan_plot(gwas_result = dts_ames2goodman_2pcs, 
                                         sig_threshold_1 = sig_threshold_1 , sig_threshold_2 = sig_threshold_2,
                                         ylim = 25,title = "DTS ~ SNP + 2 PCs (ames to 282)")



# ----------------------
# Goodman to goodman method
# ----------------------

# ------------------------------------------
#     Calculate goodman2goodman global PCs
#     uses ames2any_matrix functions
# ------------------------------------------

X_all_goodman <- combine_genofile_matrix(genofile_ames = genofile_goodman, genofile_any = genofile_goodman)

# Calculate goodman global PCs, get loadings
global_goodman_BV_3gPCs <- get_ames_BV(X_ames = genofile_goodman, Q = NULL, num_PCs = 5)

# Transfer loadings and coefficients over to all X (goodman, goodman) to get adjusted global PCs
goodman2goodman_gPCs_3gPCs_allNAM <- ames2any(X_all = X_all_goodman$X_all,
                                              Q = NULL, B = global_goodman_BV_3gPCs$B, V = global_goodman_BV_3gPCs$V)

# Subset file in half (i.e. the goodman taxa are repeated twice because that's how the function is setup)
goodman2goodman_gPCs_3gPCs_allNAM <- goodman2goodman_gPCs_3gPCs_allNAM[1:277,]
write.csv(goodman2goodman_gPCs_3gPCs_allNAM, "Q_all_3gPCs_allGoodman_goodman2goodman.csv")

# Get variance explained by first 3 pcs
# Function to get goodman B and V matrices given original X matrix
# temp <- get_ames_BV_all_pcs(genofile_goodman, num_PCs = 5)

# Load in pre-analyzed PCs
goodman2goodman_3gPCs_df <- read.csv("Q_all_3gPCs_allGoodman_goodman2goodman.csv", header = TRUE) %>% 
  select(-PC4, -PC5)
colnames(goodman2goodman_3gPCs_df) <- c("Taxa", "PC1", "PC2", "PC3")

# Map traits
map_trait(pc_df = goodman2goodman_3gPCs_df, num_covars = 3, 
          id_paste = "_fast_assoc_results_goodman2goodman_flowering_time_buckler_282_3pcs")


# Load in results (goodman2goodman)
my_files <- list.files(pattern = "_fast_assoc_results_goodman2goodman_flowering_time_buckler_282_3pcs.txt$")
my_files
all_dts <- c()
all_dts <- lapply(my_files, data.table::fread) %>% data.table::rbindlist()

# Process results
dts_goodman2goodman_3pcs <- data.frame(all_dts) %>% filter(Trait == "DTS_BLUP.Buckler_2009")

# Get significance thresholds
sig_threshold_1 = -log10(quantile(dts_goodman2goodman_3pcs$p, probs = 0.05))
sig_threshold_2 = -log10(quantile(dts_goodman2goodman_3pcs$p, probs = 0.01))

goodman2goodman_gwas_3pcs <- manhattan_plot(gwas_result = dts_goodman2goodman_3pcs, 
                                            sig_threshold_1 = sig_threshold_1, 
                                            sig_threshold_2 = sig_threshold_2,
                                            ylim = 25,
                                            title = "DTS ~ SNP + 3 PCs (282 to 282)")



# ---------------------------
# Goodman 2 goodman with 2PCs
# ---------------------------

goodman2goodman_2gPCs_df <- goodman2goodman_3gPCs_df %>% select(-PC3)

# Map traits
map_trait(pc_df = goodman2goodman_2gPCs_df, num_covars = 2, 
          id_paste = "_fast_assoc_results_goodman2goodmanflowering_time_buckler_282_2pcs")


# Load in results (goodman2goodman)
my_files <- list.files(pattern = "_fast_assoc_results_goodman2goodmanflowering_time_buckler_282_2pcs.txt$")
my_files
all_dts <- c()
all_dts <- lapply(my_files, data.table::fread) %>% data.table::rbindlist()

# Process results
dts_goodman2goodman_2pcs <- data.frame(all_dts) %>% filter(Trait == "DTS_BLUP.Buckler_2009")

# Get significance thresholds
sig_threshold_1 = -log10(quantile(dts_goodman2goodman_2pcs$p, probs = 0.05))
sig_threshold_2 = -log10(quantile(dts_goodman2goodman_2pcs$p, probs = 0.01))

goodman2goodman_gwas_2pcs <- manhattan_plot(gwas_result = dts_goodman2goodman_2pcs, 
                                            sig_threshold_1 = sig_threshold_1, 
                                            sig_threshold_2 = sig_threshold_2,
                                            ylim = 25,
                                            title = "DTS ~ SNP + 2 PCs (282 to 282)")




# Put images together
library(patchwork)
ames2goodman_gwas_2pcs | ames2goodman_gwas_3pcs | goodman2goodman_gwas_2pcs | goodman2goodman_gwas_3pcs
ggsave("two_three_pc_methods_goodman.png")



# ------------------------
# Plot by chromosome
# ------------------------

# Ames to goodman, 2 PCs
a1  <- manhattan_without_windows8(results = dts_ames2goodman_2pcs, upper_bound = 25,
                                  title = "Chrom 8: DTS ~ SNP + 2 PCs (Ames to 282)")
a2 <- manhattan_without_windows10(results = dts_ames2goodman_2pcs, upper_bound = 25,
                               title = "Chrom 10: DTS ~ SNP + 2 PCs (Ames to 282)") 

# Ames to goodman, 3 PCs
b1 <- manhattan_without_windows8(results = dts_ames2goodman_3pcs, upper_bound = 25,
                              title = "Chrom 8: DTS ~ SNP + 3 PCs (Ames to 282)")
b2 <- manhattan_without_windows10(results = dts_ames2goodman_3pcs, upper_bound = 25,
                               title = "Chrom 10: DTS ~ SNP + 3 PCs (Ames to 282)")

# Goodman to goodman, 2 PCs
c1 <- manhattan_without_windows8(results = dts_goodman2goodman_2pcs, upper_bound = 25,
                              title = "Chrom 8: DTS ~ SNP + 2 PCs (282 to 282)")
c2 <- manhattan_without_windows10(results = dts_goodman2goodman_2pcs, upper_bound = 25,
                               title = "Chrom 10: DTS ~ SNP + 2 PCs (282 to 282)")

# Goodman to goodman, 3 PCs
d1 <- manhattan_without_windows8(results = dts_goodman2goodman_3pcs, upper_bound = 25,
                              title = "Chrom 8: DTS ~ SNP + 3 PCs (282 to 282)")
d2 <- manhattan_without_windows10(results = dts_goodman2goodman_3pcs, upper_bound = 25,
                               title = "Chrom 10: DTS ~ SNP + 3 PCs (282 to 282)")

# Plot together
(a1|a2)/(b1|b2)/(c1|c2)/(d1|d2)
ggsave("chrom_8_10_two_three_pc_methods_goodman.png")



# Investigate ZNC8 peak SNP
df <- dts_ames2goodman_2pcs
df <- dts_ames2goodman_3pcs
df <- dts_goodman2goodman_2pcs
df <- dts_goodman2goodman_3pcs

sig_threshold_2 = quantile(df$p, probs = 0.01)
temp <- df %>% filter(Chr == 8 & p <= sig_threshold_2 & Pos >= 126880531-1000 & Pos <= 126880531+3000000)
max(-log10(temp$p))
temp[which.max(-log10(temp$p)),]


lala <- c(128863816,128863816 , 127830463, 127830463)



# ---------------------------
# Investigate vgt1 region
# ---------------------------

g2g_3pcs <- read.table("/workdir/mbb262/goodman282/filtered_for_pca/goodman2goodman_vgt1.txt", header = T)
a2g_1pcs <- read.table("/workdir/mbb262/goodman282/filtered_for_pca/ames2goodman_vgt1_1pc.txt", header = T)
a2g_2pcs <- read.table("/workdir/mbb262/goodman282/filtered_for_pca/ames2goodman_vgt1_2pc.txt", header = T)
a2g_3pcs <- read.table("/workdir/mbb262/goodman282/filtered_for_pca/ames2goodman_vgt1_3pc.txt", header = T)

# Variables to help plot
sig_threshold_1_1 = -log10(quantile(g2g_3pcs$p, probs = 0.05))
sig_threshold_2_1 = -log10(quantile(g2g_3pcs$p, probs = 0.01))
upper_bound_1 = 25
text_pos_1 = upper_bound-4

# Function to plot things
library(ggplot2)
library(patchwork)
man_chr8_vgt1 <- function(temp, sig_threshold_1 = sig_threshold_1_1, sig_threshold_2 = sig_threshold_2_1,
                          upper_bound = upper_bound_1, title, text_pos = text_pos_1) {
  temp <- temp %>% filter(Trait == 'DTA_BLUP;Buckler_2009')
  
  ggplot(temp, aes(x = Pos, y = -log10(p), size = -log10(p))) +
    geom_point(size = 1) +
    geom_hline(yintercept = sig_threshold_1, color = "grey40", linetype = "dashed") +
    geom_hline(yintercept = sig_threshold_2, color = "black", linetype = "dashed") +
    scale_y_continuous(expand = c(0,0), limits = c(1.5, upper_bound)) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = NULL, y = "-log10(p)", title = title) + 
    theme_minimal() +
    theme( 
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(size = 9, vjust = 0.5),
      axis.text.y = element_text(size = 9, vjust = 0.5)) +
    annotate("text", x = 135946643, y = text_pos, label = "vgt1") +
    geom_vline(xintercept=135946643) # +
    # annotate("text", x =136009216 , y = text_pos-2, label = "ZmRAP2.7") +
    # geom_vline(xintercept=136009216) 
}

# Plot
g2g_3pcs_plot <- man_chr8_vgt1(temp = g2g_3pcs, title = "Chrom 8: DTS ~ SNP + 3 PCs (282 to 282)")
a2g_1pcs_plot <- man_chr8_vgt1(temp = a2g_1pcs, title = "Chrom 8: DTS ~ SNP + 1 PCs (Ames to 282)")
a2g_2pcs_plot <- man_chr8_vgt1(temp = a2g_2pcs, title = "Chrom 8: DTS ~ SNP + 2 PCs (Ames to 282)")
a2g_3pcs_plot <- man_chr8_vgt1(temp = a2g_3pcs, title = "Chrom 8: DTS ~ SNP + 3 PCs (Ames to 282)")

g2g_3pcs_plot + a2g_1pcs_plot + a2g_2pcs_plot + a2g_3pcs_plot
ggsave("pc_methods_around_vgt1.png")


# ------------------------------
# Investigate effect estimates
# ------------------------------

e1 <- read.table("~/Downloads/glm_results/GLM_3pcs_effect_estimates.txt", header = TRUE)
r1 <- read.table("~/Downloads/glm_results/GLM_3pcs_effect_estimates_gwas_results.txt", header = TRUE)

e2 <- read.table("~/Downloads/glm_results/GLM_Stats_Filtered_flowering_buckler_282pop_vgt1_interval_maf_2pcs_effects.txt", header = TRUE)
r2 <- read.table("~/Downloads/glm_results/GLM_Stats_Filtered_flowering_buckler_282pop_vgt1_interval_maf_2pcs.txt", header = TRUE)

e3 <- read.table("~/Downloads/glm_results/GLM_Stats_hmp321_282_agpv4_merged_chr8_imputed_goodman282_vgt1_interval_vgt1_interval_maf_3pcs_effects.txt", header = TRUE)
r3 <- read.table("~/Downloads/glm_results/GLM_Stats_hmp321_282_agpv4_merged_chr8_imputed_goodman282_vgt1_interval_vgt1_interval_maf_3pcs.txt", header = TRUE)

e4 <- read.table("~/Downloads/glm_results/GLM_Stats_Filtered_flowering_buckler_282pop_vgt1_interval_maf_3pcs_nomincount_effects.txt", header = TRUE)
r4 <- read.table("~/Downloads/glm_results/GLM_Stats_Filtered_flowering_buckler_282pop_vgt1_interval_maf_3pcs_nomincount.txt", header = TRUE)

e5 <- read.table("~/Downloads/glm_results/GLM_Stats_Filtered_Filtered_Q_all_3gPCs_allGoodman + Filtered_flowering_buckler_282pop + hmp321_282_agpv4_merged_chr8_imputed_goodman282_vgt1_interval_effect.txt", header = TRUE)
r5 <- read.table("~/Downloads/glm_results/GLM_Stats_Filtered_Filtered_Q_all_3gPCs_allGoodman + Filtered_flowering_buckler_282pop + hmp321_282_agpv4_merged_chr8_imputed_goodman282_vgt1_interval.txt", header = TRUE)

e6 <- read.table("~/Downloads/glm_results/GLM_Stats_Q_all_3gPCs_allGoodman + Filtered_flowering_buckler_282pop + hmp321_282_agpv4_merged_chr8_imputed_goodman282_vgt1_interval_vgt1_interval_maf_3pcs_mincount150_effects.txt", header = TRUE)
r6 <- read.table("~/Downloads/glm_results/GLM_Stats_Q_all_3gPCs_allGoodman + Filtered_flowering_buckler_282pop + hmp321_282_agpv4_merged_chr8_imputed_goodman282_vgt1_interval_vgt1_interval_maf_3pcs_mincount150.txt", header = TRUE)

e7 <- read.table("~/Downloads/glm_results/GLM_Stats_Filtered_Q_all_3gPCs_allGoodman + Filtered_flowering_buckler_282pop + hmp321_282_agpv4_merged_chr8_imputed_goodman282_vgt1_interval_vgt1_interval_maf_3pcs_mincount100_effects.txt", header = TRUE)
r7 <- read.table("~/Downloads/glm_results/GLM_Stats_Filtered_Q_all_3gPCs_allGoodman + Filtered_flowering_buckler_282pop + hmp321_282_agpv4_merged_chr8_imputed_goodman282_vgt1_interval_vgt1_interval_maf_3pcs_mincount100.txt", header = TRUE)

e8 <- read.table("~/Downloads/glm_results/GLM_Stats_Filtered_Filtered_Q_all_3gPCs_allGoodman + Filtered_flowering_buckler_282pop + hmp321_282_agpv4_merged_chr8_imputed_goodman282_vgt1_interval_vgt1_interval_maf_3pcs_mincount55_effect.txt", header = TRUE)
r8 <- read.table("~/Downloads/glm_results/GLM_Stats_Filtered_Filtered_Q_all_3gPCs_allGoodman + Filtered_flowering_buckler_282pop + hmp321_282_agpv4_merged_chr8_imputed_goodman282_vgt1_interval_vgt1_interval_maf_3pcs_mincount55.txt", header = TRUE)

e9 <- read.table("~/Downloads/glm_results/GLM_Stats_Filtered_Q_all_3gPCs_allGoodman + Filtered_flowering_buckler_282pop + hmp321_282_agpv4_merged_chr8_imputed_goodman282_vgt1_interval_vgt1_interval_maf_3pcs_mincount55_effect.txt", header = TRUE)
r9 <- read.table("~/Downloads/glm_results/GLM_Stats_Filtered_Q_all_3gPCs_allGoodman + Filtered_flowering_buckler_282pop + hmp321_282_agpv4_merged_chr8_imputed_goodman282_vgt1_interval_vgt1_interval_maf_3pcs_mincount55.txt", header = TRUE)

e10 <- read.table("~/Downloads/glm_results/GLM_Stats_Filtered_Filtered_Q_all_3gPCs_allGoodman + Filtered_flowering_buckler_282pop + hmp321_282_agpv4_merged_chr8_imputed_goodman282_vgt1_interval_vgt1_interval_maf_124pcs_mincount55_effect.txt", header = TRUE)
r10 <- read.table("~/Downloads/glm_results/GLM_Stats_Filtered_Filtered_Q_all_3gPCs_allGoodman + Filtered_flowering_buckler_282pop + hmp321_282_agpv4_merged_chr8_imputed_goodman282_vgt1_interval_vgt1_interval_maf_124pcs_mincount55.txt", header = TRUE)

# Variables to help plot
sig_threshold_1 = -log10(quantile(r1$p, probs = 0.05))
sig_threshold_2 = -log10(quantile(r1$p, probs = 0.01))
upper_bound = 10
text_pos = upper_bound-0.5

library(dplyr)
plot_effects <- function(results, effect, title) {
  
  results <- results %>% filter(Pos >= 135731499 & Pos <= 135951643)
  effect <- effect %>% filter(Pos >= 135731499 & Pos <= 135951643)
  
  plot_results <- ggplot(results, aes(x = Pos, y = -log10(p))) +
    geom_point(size = 1) +
    geom_hline(yintercept = sig_threshold_1, color = "grey40", linetype = "dashed") +
    geom_hline(yintercept = sig_threshold_2, color = "black", linetype = "dashed") +
    scale_y_continuous(expand = c(0,0), limits = c(1.5, upper_bound)) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = NULL, y = "-log10(p)", title = paste("GWAS Results - ", title)) + 
    theme_minimal() +
    theme( 
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(size = 9, vjust = 0.5),
      axis.text.y = element_text(size = 9, vjust = 0.5)) +
    annotate("text", x = 135946643, y = text_pos, label = "vgt1") +
    geom_vline(xintercept=135946643) +
    annotate("text", x = 135736499, y = text_pos, label = "??") +
    geom_vline(xintercept=135736499)
  
  plot_effect <- ggplot(effect, aes(x = Pos, y = Estimate)) +
    geom_point(size = 1) +
    scale_y_continuous(expand = c(0,0), limits = c(-40, 40)) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = NULL, y = "Effect Estimate") + 
    theme_minimal() +
    theme( 
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(size = 9, vjust = 0.5),
      axis.text.y = element_text(size = 9, vjust = 0.5)) +
    geom_vline(xintercept=135946643) +
    geom_vline(xintercept=135736499)
  
  plot_addeffect <- ggplot(results, aes(x = Pos, y = addEffect)) +
    geom_point(size = 1) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 20)) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = NULL, y = "Additive Effect") + 
    theme_minimal() +
    theme( 
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(size = 9, vjust = 0.5),
      axis.text.y = element_text(size = 9, vjust = 0.5)) +
    geom_vline(xintercept=135946643) +
    geom_vline(xintercept=135736499)
  
  return(plot_results/plot_effect/plot_addeffect)
}


a <- plot_effects(r1, e1, "\n DTS ~ 3 PCs (Ames to 282) \n no MAF or count filter")
b <- plot_effects(r2, e2, "\n DTS ~ 2 PCs (Ames to 282) \n MAF > 0.01 & min count > 55")
c <- plot_effects(r3, e3, "\n DTS ~ 3 PCs (Ames to 282) \n MAF > 0.01 & min count > 55")
d <- plot_effects(r4, e4, "\n DTS ~ 3 PCs (Ames to 282) \n MAF > 0.01 & no min count filter")
e <- plot_effects(r5, e5, "\n DTS ~ 2 PCs (Ames to 282) \n no MAF or count filter")
f <- plot_effects(r6, e6, "\n DTS ~ 3 PCs (Ames to 282) \n MAF > 0.01 & min count > 150")
g <- plot_effects(r7, e7, "\n DTS ~ 3 PCs (Ames to 282) \n MAF > 0.01 & min count > 55 & class 5")
h <- plot_effects(r8, e8, "\n DTS ~ 2 PCs (Ames to 282) \n MAF > 0.01 & min count > 55 & class 5")
i <- plot_effects(r9, e9, "\n DTS ~ PCs 1,2,4,5 (Ames to 282) \n MAF > 0.01 & min count > 55 & class 5")
j <- plot_effects(r10, e10, "\n DTS ~ PCs 1,2, 4 (Ames to 282) \n MAF > 0.01 & min count > 55 & class 5")

g|i|j
a|e|b|c|d|f

b | c
ggsave("effects.png")


focal_2pcs = 135736499
focal_3pcs = 135942722

start = 135731499
end = start+5000
temp <- r3 %>% filter(Pos >= start & Pos <= end)
temp_r2 <- r2 %>% filter(Pos >= start & Pos <= end)

r7[844,]
r8 %>% filter(Pos == 135736499)
