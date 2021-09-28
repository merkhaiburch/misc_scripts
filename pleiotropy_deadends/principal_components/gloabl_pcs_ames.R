# ---------------------------------------
# Merritt Burch
# mbb262@cornell.edu
# 2019-09-09
# Script to do PCA and export PCs
#
# Calculate PCs on mean imputed data
# (vcf, distance matrix and kinship matrix)
# ---------------------------------------

# Set directory
setwd("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/ames_tests/Sep9_Ames_imputationTest")

# read in imputed data
impute_vcf <- read.table("impute_output1.txt", header = TRUE)
impute_distance <- read.table("ZeaGBSv27_20171204_AGPv4_Ames_numerical_imputation_genotypeDistance.txt", header = FALSE)
impute_kinship <- read.table("ZeaGBSv27_20171204_AGPv4_Ames_numerical_imputation_kinshipMatrix.txt", header = FALSE)


# ----------------
# Calculate PCs
# ----------------

# Calculate PCs (imputed)
impute_distance_pca <- prcomp(na.omit(impute_distance[,2:ncol(impute_distance)]), scale = FALSE)
impute_kinship_pca <- prcomp(na.omit(impute_kinship[,2:ncol(impute_kinship)]), scale = FALSE)


# ----------------
# Extract out PCs
# ----------------

# Extract out PCs (imputed)
impute_distance_loadings <- as.data.frame(impute_distance_pca$rotation)
impute_kinship_loadings <- as.data.frame(impute_kinship_pca$rotation)


# --------------------------------
# Run models
# Run with all 17 Ames phenotypes
# --------------------------------

# Add taxa names to PCA loading score data.frames
impute_distance_loadings <- cbind(impute_distance$Taxa, impute_distance_loadings)
impute_kinship_loadings <- cbind(impute_kinship$Taxa, impute_kinship_loadings)

# Fix column name in new data.frame
colnames(impute_distance_loadings)[1] <- "Taxa"
colnames(impute_kinship_loadings)[1] <- "Taxa"

# Load all ames phenotypes and subset down to only what I want
ames_phenos <- read.table("/Users/mbb262-admin/git_projects/haplotype_v_snp_gwas/data/all_Ames_Phenos.txt", header = TRUE)
ames_phenos <- ames_phenos[,c(1,6:22)]

# Format names for merging in imputed files
impute_vcf$simpleNames <- gsub(":[0-9]{9}", "", tolower(impute_vcf$Taxa))
impute_distance_loadings$simpleNames <- gsub(":[0-9]{9}", "", tolower(impute_distance_loadings$Taxa))
impute_kinship_loadings$simpleNames <- gsub(":[0-9]{9}", "", tolower(impute_kinship_loadings$Taxa))

# Merge phenotypes with PCs
impute_vcf_mergedPhenos <- merge(y = impute_vcf, x = ames_phenos, by.y = "simpleNames", by.x = "Entry_ID")
impute_distance_mergedPhenos <- merge(y = impute_distance_loadings, x = ames_phenos, by.y = "simpleNames", by.x = "Entry_ID")
impute_kinship_mergedPhenos <- merge(y = impute_kinship_loadings, x = ames_phenos, by.y = "simpleNames", by.x = "Entry_ID")


# ------------------------------------------------------
# Run model
# y (17 ames phenotypes) ~ 5 PCs (vcf/distance/kinship)
# ------------------------------------------------------

# Make a loop that collects model p-values for all 17 traits
library(broom)
library(magrittr)

pvals <- data.frame()
all_anovas <- function(model_df){
  for (i in seq(2,18)){
    model <- lm(model_df[,i] ~ PC1+PC2+PC3+PC4+PC5, model_df)
    tmp <- anova(model) %>% broom::tidy()
    pvals <- rbind(pvals, t(tmp$p.value))
  }
  traits <- data.frame(colnames(impute_vcf_mergedPhenos[,2:18]))
  pvals_n_traits <- cbind(traits, pvals[,-6])
  colnames(pvals_n_traits) <- c("trait", "PC1", "PC2", "PC3", "PC4", "PC5")
  return(pvals_n_traits)
}


# Run the loops
vcf_lm_pvals <- all_anovas(impute_vcf_mergedPhenos)
distance_lm_pvals <- all_anovas(impute_distance_mergedPhenos)
kinship_lm_pvals <- all_anovas(impute_kinship_mergedPhenos)


# -------------
# Plot results
# -------------

# On vcfs
boxplot(vcf_lm_pvals[,-1],
        xlab = "First 5 PCs",
        ylab = "Model p-value",
        main = "y(17 Ames phenotypes) ~ first 5 PCs(mean imputation on SNPs directly)")

# on distance matrix
boxplot(distance_lm_pvals[,-1],
        xlab = "First 5 PCs",
        ylab = "Model p-value",
        main = "y(17 Ames phenotypes) ~ first 5 PCs(mean imputation on distance matrix)")

# On kinship matrix
boxplot(kinship_lm_pvals[,-1],
        xlab = "First 5 PCs",
        ylab = "Model p-value",
        main = "y(17 Ames phenotypes) ~ first 5 PCs(mean imputation on kinship matrix)")


# -------------------
# Plot cumulative sum
# -------------------

par(mfrow = c(1,3))

vcf_var <- read.table("impute_output2.txt", header = TRUE)
plot(vcf_var$cumulativeproportion,
     main = "Cumulative variance explained (PCs on vcf)")

distance_var <- cumsum(impute_distance_pca$sdev^2/sum(impute_distance_pca$sdev^2))
plot(distance_var,
     main = "Cumulative variance explained (PCs on distance matrix)")

kinship_var <- cumsum(impute_kinship_pca$sdev^2/sum(impute_kinship_pca$sdev^2))
plot(kinship_var,
     main = "Cumulative variance explained (PCs on kinship matrix)")


# Make nice table

# Make into dataframe
vcf_var <- data.frame(vcf_var)
distance_var <- data.frame(distance_var)
kinship_var <- data.frame(kinship_var)

# Add in IDs
vcf_var$Method <- rep("PCs on VCF", nrow(vcf_var))
distance_var$Method <- rep("PCs on Distance", nrow(distance_var))
kinship_var$Method <- rep("PCs on Kinship", nrow(kinship_var))

# Add variance explained by each PC
distance_var$proportionoftotal <- impute_distance_pca$sdev^2/sum(impute_distance_pca$sdev^2)
kinship_var$proportionoftotal <- impute_kinship_pca$sdev^2/sum(impute_kinship_pca$sdev^2)

# Subset
vcf_var_subset <- vcf_var[1:5,3:5]
distance_var_subset <- distance_var[1:5,c(3,1,2)]
kinship_var_subset <- kinship_var[1:5, c(3,1,2)]

# Add in PC names
vcf_var_subset <- cbind(c("PC1", "PC2", "PC3", "PC4", "PC5") , vcf_var_subset)
distance_var_subset <- cbind(c("PC1", "PC2", "PC3", "PC4", "PC5") , distance_var_subset)
kinship_var_subset <- cbind(c("PC1", "PC2", "PC3", "PC4", "PC5") , kinship_var_subset)

# Rename columns
colnames(vcf_var_subset) <- c("PC", "proportionoftotal", "cumulativeproportion", "Method")
colnames(distance_var_subset) <- c("PC", "proportionoftotal", "cumulativeproportion", "Method")
colnames(kinship_var_subset) <- c("PC", "proportionoftotal", "cumulativeproportion", "Method")

# combine
together <- rbind(vcf_var_subset, distance_var_subset, kinship_var_subset)
temp <- as.data.frame(t(together))
write.csv(together, "top_5PCs_variance_mean_imputation_vcf_distance_kinship.csv",
            row.names = FALSE, quote = FALSE)


# ----------------------------
# Plot PC axes with annotation
# ----------------------------

# Load packages
library(ggplot2)

# load supp. figure 1 from Romay 2013
romay <- read.csv("/Users/mbb262-admin/Box Sync/Cornell_PhD/labProjects/hap_gwa/ames_tests/Romay_2013_supp_table1.csv", header = TRUE)
romay$Accesion.N <- tolower(romay$Accesion.N)

# Merge pca loadings with romay annotations
vcf_loadings_romay <- merge(x = impute_vcf_mergedPhenos, by.x = "Entry_ID",
                            y = romay, by.y = "Accesion.N")
distance_loadings_romay <- merge(x = impute_distance_mergedPhenos, by.x = "Entry_ID",
                                 y = romay, by.y = "Accesion.N")
kinship_loadings_romay <- merge(x = impute_kinship_mergedPhenos, by.x = "Entry_ID",
                                 y = romay, by.y = "Accesion.N")



# PC1 by PC2
ggplot(vcf_loadings_romay, aes(x=PC1, y=PC2, color = Pop.structure)) + 
  geom_point() +
  ggtitle("Plot of PCs on mean imputed vcf (PCs directly on SNPs)")

ggplot(distance_loadings_romay, aes(x=PC1, y=PC2, color = Pop.structure)) + 
  geom_point() +
  ggtitle("Plot of PCs on mean imputed distance matrix")

ggplot(kinship_loadings_romay, aes(x=PC1, y=PC2, color = Pop.structure)) + 
  geom_point() +
  ggtitle("Plot of PCs on mean imputed kiship matrix")


# PC2 by PC3
ggplot(vcf_loadings_romay, aes(x=PC2, y=PC3, color = Pop.structure)) + 
  geom_point() +
  ggtitle("Plot of PCs on mean imputed vcf (PCs directly on SNPs)")

ggplot(distance_loadings_romay, aes(x=PC2, y=PC3, color = Pop.structure)) + 
  geom_point() +
  ggtitle("Plot of PCs on mean imputed distance matrix")

ggplot(kinship_loadings_romay, aes(x=PC2, y=PC3, color = Pop.structure)) + 
  geom_point() +
  ggtitle("Plot of PCs on mean imputed kiship matrix")


# PC3 by PC4
ggplot(vcf_loadings_romay, aes(x=PC3, y=PC4, color = Pop.structure)) + 
  geom_point() +
  ggtitle("Plot of PCs on mean imputed vcf (PCs directly on SNPs)")

ggplot(distance_loadings_romay, aes(x=PC3, y=PC4, color = Pop.structure)) + 
  geom_point() +
  ggtitle("Plot of PCs on mean imputed distance matrix")

ggplot(kinship_loadings_romay, aes(x=PC3, y=PC4, color = Pop.structure)) + 
  geom_point() +
  ggtitle("Plot of PCs on mean imputed kiship matrix")


# PC4 by PC5
ggplot(vcf_loadings_romay, aes(x=PC4, y=PC5, color = Pop.structure)) + 
  geom_point() +
  ggtitle("Plot of PCs on mean imputed vcf (PCs directly on SNPs)")

ggplot(distance_loadings_romay, aes(x=PC4, y=PC5, color = Pop.structure)) + 
  geom_point() +
  ggtitle("Plot of PCs on mean imputed distance matrix")

ggplot(kinship_loadings_romay, aes(x=PC4, y=PC5, color = Pop.structure)) + 
  geom_point() +
  ggtitle("Plot of PCs on mean imputed kiship matrix")
