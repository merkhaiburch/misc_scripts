# ---------------------------------------
# Merritt Burch
# mbb262@cornell.edu
# 2019-09-19
# Script to do PCA, plot them, and do
# analyses
#
# Calculate PCs on imputed data --> (mean, beagle, kinship)
# Look at PC calculation inputs --> (vcf, distance matrix and kinship matrix)
# Data filtering: table sites: MinCount = 2500, MAF 0.0001
# table taxa: MinFreq = 0.55
# ---------------------------------------

# Set directory
setwd("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/ames_tests/Sep20_whole_pipeline")

# -------------------------
# Load in all data at once
# -------------------------

# Load all ames phenotypes and subset down to only what I want
ames_phenos <- read.table("/Users/mbb262-admin/git_projects/haplotype_v_snp_gwas/data/all_Ames_Phenos.txt", header = TRUE)
ames_phenos <- ames_phenos[,c(1,6:22)]

# No imputation
no_imputation_vcf_pcs <- read.table("no_imputation_output1.txt", header = TRUE)
no_imputation_vcf_eigen <- read.table("no_imputation_output2.txt",  header = TRUE)
no_imputation_distance <- read.table("ZeaGBSv27_20171204_AGPv4_Ames_no_imputation_genotypeDistance.txt", header = FALSE)
no_imputation_kinship <- read.table("ZeaGBSv27_20171204_AGPv4_Ames_no_imputation_kinshipMatrix.txt", header = FALSE)

# Mean imputation
mean_vcf_pcs <- read.table("mean_impute_output1.txt", header = TRUE)
mean_vcf_eigen <- read.table("mean_impute_output2.txt",  header = TRUE)

# Beagle imputation
beagle_vcf_pcs <- read.table("beagle_impute_output1.txt", header = TRUE)
beagle_vcf_eigen <- read.table("beagle_impute_output2.txt", header = TRUE)
beagle_distance <- read.table("ZeaGBSv27_20171204_AGPv4_Ames_beagle_imputation_genotypeDistance.txt", header = FALSE)
beagle_kinship <- read.table("ZeaGBSv27_20171204_AGPv4_Ames_beagle_imputation_kinshipMatrix.txt", header = FALSE)

# KNNI imputation
knni_vcf_pcs <- read.table("knni_impute_output1.txt", header = TRUE)
knni_vcf_eigen <- read.table("knni_impute_output2.txt", header = TRUE)
knni_distance <- read.table("ZeaGBSv27_20171204_AGPv4_Ames_knni_imputation_genotypeDistance.txt", header = FALSE)
knni_kinship <- read.table("ZeaGBSv27_20171204_AGPv4_Ames_knni_imputation_kinshipMatrix.txt", header = FALSE)

# fillin imputation
fillin_vcf_pcs <- read.table("fillin_impute_output1.txt", header = TRUE)
fillin_vcf_eigen <- read.table("fillin_impute_output2.txt", header = TRUE)
fillin_distance <- read.table("ZeaGBSv27_20171204_AGPv4_Ames_fillin_imputation_genotypeDistance.txt", header = FALSE)
fillin_kinship <- read.table("ZeaGBSv27_20171204_AGPv4_Ames_fillin_imputation_kinshipMatrix.txt", header = FALSE)


# ---------------------------------------
# Do PCA on distance and kinship matrices
# ---------------------------------------

# Use this one it provides the right output
decomp <- function(G, k=2, center=TRUE) {
  
  # stopifnot(isSymmetric(G))
  
  n <- nrow(G)
  
  # Centering (in case it is not already)
  if (center) {
    H <- diag(n) - matrix(1, nrow=n, ncol=n)/n
    G <- H %*% G %*% H
  }
  
  # Decomposition
  ED <- eigen(G)
  
  # Output PCs
  return(list(PC=ED$vectors[, 1:k] %*% diag(sqrt(ED$values[1:k])),
              variance=ED$values[1:k],
              variance_ratio=ED$values[1:k]/sum(ED$values)))
  
}

# No imputation
# imputation on vcf file done in tassel
no_imputation_distance_pcs <- cmdscale(as.dist(no_imputation_distance[,-1]), k=5, eig = T)
no_imputation_kinship_pcs <- decomp(as.matrix(no_imputation_kinship[,-1]), k=5)

# Mean imputation
# Nothing to do here folks

# Beagle imputation
# imputation on vcf file done in tassel
beagle_distance_pcs <- cmdscale(as.dist(beagle_distance[,-1]), k=5, eig = T)
beagle_kinship_pcs <- decomp(as.matrix(beagle_kinship[,-1]), k=5)

# KNNI imputation
# imputation on vcf file done in tassel
knni_distance_pcs <- cmdscale(as.dist(knni_distance[,-1]), k=5, eig = T)
knni_kinship_pcs <- decomp(as.matrix(knni_kinship[,-1]), k=5)

# FILLIN imputation
# imputation on vcf file done in tassel
fillin_distance_pcs <- cmdscale(as.dist(fillin_distance[,-1]), k=5, eig = T)
fillin_kinship_pcs <- decomp(as.matrix(fillin_kinship[,-1]), k=5)


# ------------------------------------------
# Get variance explained by each PC (table)
# ------------------------------------------

# No imputation
# variance by vcf already done
no_imputation_distance_pcs_var <- 100 * no_imputation_distance_pcs$eig / sum(no_imputation_distance_pcs$eig)
no_imputation_kinship_pcs_var <- data.frame(no_imputation_kinship_pcs$variance_ratio*100)

# Mean Imputation
# Already loaded in 'mean_vcf_eigen'

# Beagle
# variance by vcf already done
beagle_distance_pcs_var <- 100 * beagle_distance_pcs$eig / sum(beagle_distance_pcs$eig)
beagle_kinship_pcs_var <- data.frame(beagle_kinship_pcs$variance_ratio*100)

# KNNI
# variance by vcf already done
knni_distance_pcs_var <- 100 * knni_distance_pcs$eig / sum(knni_distance_pcs$eig)
knni_kinship_pcs_var <- data.frame(knni_kinship_pcs$variance_ratio*100)

#FILLIN
# variance by vcf already done
fillin_distance_pcs_var <- 100 * fillin_distance_pcs$eig / sum(fillin_distance_pcs$eig)
fillin_kinship_pcs_var <- data.frame(fillin_kinship_pcs$variance_ratio*100)


# Functions to format names and PCs and add them to the data.frame
comb_vcf <- function(df_vcf, id_string, imp_string){
  temp <- data.frame(df_vcf[1:5,3]*100)
  temp$pc <-seq(1,5)
  temp$method <- rep(id_string, 5)
  temp$imputation <- rep(imp_string,5)
  temp$input <- rep("vcf",5)
  colnames(temp) <- c("proportion_variance", "pc", "method", "imputation", "input")
  temp <- temp[,c(2,1,3,4,5)]
  return(temp)
}

comb_dist <- function(df_dist, id_string, imp_string){
  temp <- data.frame(df_dist[1:5])
  temp$pc <-seq(1,5)
  temp$method <- rep(id_string, 5)
  temp$imputation <- rep(imp_string,5)
  temp$input <- rep("distance",5)
  colnames(temp) <- c("proportion_variance", "pc", "method", "imputation", "input")
  temp <- temp[,c(2,1,3,4,5)]
  return(temp)
}

comb_kin <- function(df_kin, id_string, imp_string){
  temp <- data.frame(df_kin)
  temp$pc <-seq(1,5)
  temp$method <- rep(id_string, 5)
  temp$imputation <- rep(imp_string,5)
  temp$input <- rep("kinship",5)
  colnames(temp) <- c("proportion_variance", "pc", "method", "imputation", "input")
  temp <- temp[,c(2,1,3,4,5)]
  return(temp)
}

# Make the combined table
first_5 <- rbind(comb_vcf(no_imputation_vcf_eigen, "no imputation - directly on vcf", "no imputation"),
             comb_dist(no_imputation_distance_pcs_var, "no imputation - distance matrix", "no imputation"),
             comb_kin(no_imputation_kinship_pcs_var, "no imputation - kinship matrix", "no imputation"),
             
             comb_vcf(mean_vcf_eigen, "mean imputation - directly on vcf", "mean imputation"),
             
             comb_vcf(beagle_vcf_eigen, "beagle - directly on vcf", "beagle"),
             comb_dist(beagle_distance_pcs_var, "beagle - distance matrix", "beagle"),
             comb_kin(beagle_kinship_pcs_var, "beagle - kinship matrix", "beagle"),
             
             comb_vcf(knni_vcf_eigen, "knni - directly on vcf", "KNNi"),
             comb_dist(knni_distance_pcs_var, "knni - distance matrix", "KNNi"),
             comb_kin(knni_kinship_pcs_var, "knni - kinship matrix", "KNNi"),
             
             comb_vcf(fillin_vcf_eigen, "fillin - directly on vcf", "fillin"),
             comb_dist(fillin_distance_pcs_var, "fillin - distance matrix", "fillin"),
             comb_kin(fillin_kinship_pcs_var, "fillin - kinship matrix", "fillin"))

# format dataframe
first_5$imputation <- as.factor(first_5$imputation)
first_5$input  <- as.factor(first_5$input)
first_5$pc <- as.numeric(first_5$pc)

# Plot it out
library(ggplot2)
ggplot(first_5, aes(x=pc, y=proportion_variance, colour=input, shape = imputation)) + 
  geom_point(size=3) + 
  geom_line() +
  xlab("PC") + 
  ylab("Percent Variance Explained")


# --------------------------------------
# Correlation matrix between forst 5 PCs
# --------------------------------------

# Subset out first 5 PCs from each method,
# merge them into one dataframe, 
# add useful names
# Plot

# Note: No imputation, mean/numeric imputation, knni, and fillin all had same imput file and outputted analysis
# on file that was n=3465
# Beagle had the same exact input but ended up being n=3460. I don't understand why

# beagle_vcf_pcs[,2:6],
# data.frame(beagle_distance_pcs$points[,1:5]), 
# data.frame(beagle_kinship_pcs[,1:5])

# data.frame to collect all PCs
first_5_pcs_all <- cbind(no_imputation_vcf_pcs[,-c(1,7)],
                         data.frame(no_imputation_distance_pcs$points[,1:5]),
                         data.frame(no_imputation_kinship_pcs$PC[,1:5]),
                                       
                         mean_vcf_pcs[,-c(1,7)],
                                       
                         knni_vcf_pcs[,2:6], 
                         data.frame(knni_distance_pcs$points[,1:5]),
                         data.frame(knni_kinship_pcs$PC[,1:5]),
                         
                         fillin_vcf_pcs[,2:6], 
                         data.frame(fillin_distance_pcs$points[,1:5]),
                         data.frame(fillin_kinship_pcs$PC[,1:5]))

colnames(first_5_pcs_all) <- c("no_imputation_vcf_pc1","no_imputation_vcf_pc2","no_imputation_vcf_pc3","no_imputation_vcf_pc4","no_imputation_vcf_pc5",
                                "no_imputation_dist_pc1","no_imputation_dist_pc2","no_imputation_dist_pc3","no_imputation_dist_pc4","no_imputation_dist_pc5",
                                "no_imputation_kin_pc1","no_imputation_kin_pc2","no_imputation_kin_pc3","no_imputation_kin_pc4","no_imputation_kin_pc5",
                                
                                "mean_vcf_pc1","mean_vcf_pc2","mean_vcf_pc3","mean_vcf_pc4","mean_vcf_pc5",
                                
                                "knni_vcf_pc1","knni_vcf_pc2","knni_vcf_pc3","knni_vcf_pc4","knni_vcf_pc5",
                                "knni_dist_pc1","knni_dist_pc2","knni_dist_pc3","knni_dist_pc4","knni_dist_pc5",
                                "knni_kin_pc1","knni_kin_pc2","knni_kin_pc3","knni_kin_pc4","knni_kin_pc5",
                                
                                "fillin_vcf_pc1","fillin_vcf_pc2","fillin_vcf_pc3","fillin_vcf_pc4","fillin_vcf_pc5",
                                "fillin_dist_pc1","fillin_dist_pc2","fillin_dist_pc3","fillin_dist_pc4","fillin_dist_pc5",
                                "fillin_kin_pc1","fillin_kin_pc2","fillin_kin_pc3","fillin_kin_pc4","fillin_kin_pc5")

# Load packages
library(ggcorrplot)

# Compute a correlation matrix
corr <- round(cor(first_5_pcs_all), 1)

# Compute a matrix of correlation p-values
p.mat <- cor_pmat(corr)

# Correlation matrix of all methods
ggcorrplot(corr, hc.order = F, type = "lower",
           lab = TRUE)


# ------------------------------------------------------
# Run model --> collect R2 values, not p values
# y (17 ames phenotypes) ~ 5 PCs (vcf/distance/kinship)
# ------------------------------------------------------

# Function merges phenotypes with PCs, runs model,
# returns R^2 value in dataframe

# Load packages
library(broom)
library(magrittr)

# Run loop
r_vals <- data.frame()
all_anovas <- function(model_df_pcs, method_name, input){
  model_df_pcs$Taxa <- gsub(":[0-9]{9}", "", tolower(model_df_pcs$Taxa))
  merged_data <- merge(y = model_df_pcs, x = ames_phenos,
                       by.y = "Taxa", by.x = "Entry_ID")
  for (i in seq(2,18)){
    model <- lm(merged_data[,i] ~ PC1+PC2+PC3+PC4+PC5, merged_data)
    tmp <- data.frame(colnames(merged_data)[i], summary(model)$r.squared, summary(model)$adj.r.squared)
    tmp <- cbind(tmp, method_name, input)
    r_vals <- rbind(r_vals, tmp)
  }
  colnames(r_vals) <- c("phenotype", "r_squared", "adj_r_squared", "method", "input")
  return(r_vals)
}

# For distance matrices 
r_vals <- data.frame()
all_anovas_distance <- function(model_df_pcs, vcf_compliment, method_name, input){
  # Extract out firist 5 pcs
  get_5_pcs <- data.frame(model_df_pcs$points[,1:5])
  
  # Combine taxa IDs with PCs
  get_5_pcs$Taxa <-vcf_compliment$Taxa
  
  # Reorder columns and rename
  get_5_pcs <- get_5_pcs[,c(6,1:5)]
  colnames(get_5_pcs) <- c("Taxa", "PC1", "PC2", "PC3", "PC4", "PC5")
  
  # Format taxa names
  get_5_pcs$Taxa <- gsub(":[0-9]{9}", "", tolower(get_5_pcs$Taxa))
  merged_data <- merge(y = get_5_pcs, x = ames_phenos,
                       by.y = "Taxa", by.x = "Entry_ID")
  for (i in seq(2,18)){
    model <- lm(merged_data[,i] ~ PC1+PC2+PC3+PC4+PC5, merged_data)
    tmp <- data.frame(colnames(merged_data)[i], summary(model)$r.squared, summary(model)$adj.r.squared)
    tmp <- cbind(tmp, method_name, input)
    r_vals <- rbind(r_vals, tmp)
  }
  colnames(r_vals) <- c("phenotype", "r_squared", "adj_r_squared", "method", "input")
  return(r_vals)
}

# For kinship matrices
r_vals <- data.frame()
all_anovas_kinship <- function(model_df_pcs, vcf_compliment, method_name, input){
  # Extract out firist 5 pcs
  get_5_pcs <- data.frame(model_df_pcs$PC[,1:5])
  
  # Combine taxa IDs with PCs
  get_5_pcs$Taxa <-vcf_compliment$Taxa
  
  # Reorder columns and rename
  get_5_pcs <- get_5_pcs[,c(6,1:5)]
  colnames(get_5_pcs) <- c("Taxa", "PC1", "PC2", "PC3", "PC4", "PC5")
  
  # Format taxa names
  get_5_pcs$Taxa <- gsub(":[0-9]{9}", "", tolower(get_5_pcs$Taxa))
  merged_data <- merge(y = get_5_pcs, x = ames_phenos,
                       by.y = "Taxa", by.x = "Entry_ID")
  for (i in seq(2,18)){
    model <- lm(merged_data[,i] ~ PC1+PC2+PC3+PC4+PC5, merged_data)
    tmp <- data.frame(colnames(merged_data)[i], summary(model)$r.squared, summary(model)$adj.r.squared)
    tmp <- cbind(tmp, method_name, input)
    r_vals <- rbind(r_vals, tmp)
  }
  colnames(r_vals) <- c("phenotype", "r_squared", "adj_r_squared", "method", "input")
  return(r_vals)
}

# Run function for all vcf and distance matrices
get_r_values <- rbind(all_anovas(no_imputation_vcf_pcs, "no_imputation", "vcf"),
                      all_anovas_distance(no_imputation_distance_pcs, vcf_compliment = no_imputation_vcf_pcs, "no_imputation", "distance"),
                      all_anovas_kinship(no_imputation_kinship_pcs, vcf_compliment = no_imputation_vcf_pcs, "no_imputation", "kinship"),
                      
                      all_anovas(mean_vcf_pcs, "mean", "vcf"),
                      
                      all_anovas(beagle_vcf_pcs, "beagle", "vcf"),
                      all_anovas_distance(beagle_distance_pcs, vcf_compliment = beagle_vcf_pcs, "beagle", "distance"),
                      all_anovas_kinship(beagle_kinship_pcs, vcf_compliment = beagle_vcf_pcs, "beagle", "kinship"),
                      
                      all_anovas(knni_vcf_pcs, "knni", "vcf"),
                      all_anovas_distance(knni_distance_pcs, vcf_compliment = knni_vcf_pcs, "knni", "distance"),
                      all_anovas_kinship(knni_kinship_pcs, vcf_compliment = knni_vcf_pcs, "knni", "kinship"),
                      
                      all_anovas(fillin_vcf_pcs, "fillin", "vcf"),
                      all_anovas_distance(fillin_distance_pcs, vcf_compliment = fillin_vcf_pcs, "fillin", "distance"),
                      all_anovas_kinship(fillin_kinship_pcs, vcf_compliment = fillin_vcf_pcs, "fillin", "kinship"))

# Make a new column combinig method and input
get_r_values$method_input <- paste(get_r_values$method, get_r_values$input, sep = " - ")

# Plot put results
ggplot(get_r_values, aes(x=method_input, y=adj_r_squared, fill=method)) + 
  geom_boxplot() +
  labs(title="Plot of R^2 values- Model: y (17 ames phenotypes) ~ 5 PCs (vcf/distance/kinship)",
       x="Method", y = "Length") +
  theme(axis.text.x = element_text(angle = 75, vjust = 0.5, hjust=0.5))





