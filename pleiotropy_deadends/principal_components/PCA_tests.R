# ---------------------------------------
# Merritt Burch
# mbb262@cornell.edu
# 2019-07-24
# Script to visulize TASSEL GWAS results
# ---------------------------------------

# Set workdir
setwd("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/pcs_fastAssociation_282")

# ggman
library(ggman)


# ------------
#   25 PCs
# Plotting Fast-Association results from TASSEL
# using PCs as covariates in the 282 test data
# ------------

# Import data (have 3 iterations with 25, 50, and 100 PCs)
results <- read.table("Fast Association_Filtered_BLUEs_mdp_phenotype + PC_25_mdp_genotype + mdp_genotype.txt",
                      header = TRUE)

# Subset down the file to be compatible with qqman
dta <- results[which(results$Trait == "dpoll"), c(2:4,7)]
earHT <- results[which(results$Trait == "EarHT"), c(2:4,7)]

# 25 PCs
png("dta_282_25PCs.png",width = 5, height = 3.5, units = "in", res = 500)
ggman(dta, chrom="Chr", bp="Pos", snp="Marker", pvalue="p", pointSize = 0.3, lineColour = "black",
      sigLine = -log10(0.05), title = "Days to Anthesis in 282 (y~x+25 PCs)")
dev.off()

png("earHT_282_25PCs.png",width = 5, height = 3.5, units = "in", res = 500)
ggman(earHT, chrom="Chr", bp="Pos", snp="Marker", pvalue="p", pointSize = 0.3, lineColour = "black",
      sigLine = -log10(0.05), title = "Days to Silking in 282 (y~x+25 PCs)")
dev.off()


# ---------
#  50 PCs
# Plotting Fast-Association results from TASSEL
# using PCs as covariates in the 282 test data
# --------- 

# Import data (have 3 iterations with 25, 50, and 100 PCs)
results_50 <- read.table("Fast Association_Filtered_BLUEs_mdp_phenotype + PC_50_mdp_genotype + mdp_genotype.txt",
                      header = TRUE)

# Subset down the file to be compatible with qqman
dta_50 <- results_50[which(results_50$Trait == "dpoll"), c(2:4,7)]
earHT_50 <- results_50[which(results_50$Trait == "EarHT"), c(2:4,7)]

# 50 PCs
png("dta_282_50PCs.png",width = 5, height = 3.5, units = "in", res = 500)
ggman(dta_50, chrom="Chr", bp="Pos", snp="Marker", pvalue="p", pointSize = 0.3, lineColour = "black",
      sigLine = -log10(0.05), title = "Days to Anthesis in 282 (y~x+50 PCs)")
dev.off()

png("earHT_282_50PCs.png",width = 5, height = 3.5, units = "in", res = 500)
ggman(earHT_50, chrom="Chr", bp="Pos", snp="Marker", pvalue="p", pointSize = 0.3, lineColour = "black",
      sigLine = -log10(0.05), title = "Days to Silking in 282 (y~x+50 PCs)")
dev.off()


# --------
# 100 PCs
# Plotting Fast-Association results from TASSEL
# using PCs as covariates in the 282 test data
# --------

# Import data (have 3 iterations with 25, 50, and 100 PCs)
results_100 <- read.table("Fast Association_Filtered_BLUEs_mdp_phenotype + PC_100_mdp_genotype + mdp_genotype.txt",
                          header = TRUE)

# Subset down the file to be compatible with qqman
dta_100 <- results_100[which(results_100$Trait == "dpoll"), c(2:4,7)]
earHT_100 <- results_100[which(results_100$Trait == "EarHT"), c(2:4,7)]

# 50 PCs
png("dta_282_100PCs.png",width = 5, height = 3.5, units = "in", res = 500)
ggman(dta_100, chrom="Chr", bp="Pos", snp="Marker", pvalue="p", pointSize = 0.3, lineColour = "black",
      sigLine = -log10(0.05), title = "Days to Anthesis in 282 (y~x+100 PCs)")
dev.off()

png("earHT_282_100PCs.png",width = 5, height = 3.5, units = "in", res = 500)
ggman(earHT_100, chrom="Chr", bp="Pos", snp="Marker", pvalue="p", pointSize = 0.3, lineColour = "black",
      sigLine = -log10(0.05), title = "Days to Silking in 282 (y~x+100 PCs)")
dev.off()

N <- nrow(dta)

expected_pvals <- -log10( seq(1/N,1,by=1/N) )

plot(expected_pvals, sort(dta$p, decreasing = TRUE),
     xlab='expected-pvals (-log10)',
     ylab='observed-pvals (-log10)',
     main='QQ plot (basic)' )

abline(a=0,b=1, col='red')


# -----------------------------------------
# Plotting one nam family by top 5 PCs
# PCA on Genotypes
# `/iplant/home/shared/panzea/genotypes/GBS/v27/ZeaGBSv27_publicSamples_imputedV5_AGPv4-181023.vcf.gz`
# ----------------------------------------

# read genotype table
geno <- read.table("ZeaGBSv27_publicSamples_imputedV5_AGPv4-181023_NAM_Filter400kMarkers_FilterKy21.txt", header = T)

# Genotyoes (only numeric)
geno2 <- geno[,2:length(geno)]

# Turn genotypes into values (recheck this, think I'm setting values strangely)
geno3 <- as.data.frame(lapply(geno2, function(x) as.numeric(x)))

# Transpose markers
geno3 <- t(geno3)

# Run PCA on markers
results_pca <- prcomp(geno3, scale = FALSE)

# Extract out PCs
loadings <- as.data.frame(results_pca$rotation)

# Add in site information
siteNames <- read.table("position_list_ky21.txt", header = TRUE)
loadings <- cbind(siteNames[,1:3],loadings)

# Plot pcas by site
plot(loadings$`siteNames[, 2]`, loadings$PC1)

# ggplot, color by chromosome
library(ggplot2)

# PC1
p1 <- ggplot(loadings, aes(Site, abs(PC1))) +
   geom_point(aes(colour = factor(Chromosome))) +
   ggtitle("Plot of PC1 by genomic position in Z014 (B73xKy21)") +
   xlab("Genomic Position")

# PC2
p2 <- ggplot(loadings, aes(Site, PC2)) +
   geom_point(aes(colour = factor(Chromosome))) +
   ggtitle("Plot of PC2 by genomic position in Z014 (B73xKy21)") +
   xlab("Genomic Position")
is.na(x)
# PC3
p3 <- ggplot(loadings, aes(Site, PC3)) +
   geom_point(aes(colour = factor(Chromosome))) +
   ggtitle("Plot of PC3 by genomic position in Z014 (B73xKy21)") +
   xlab("Genomic Position")

# PC4
p4 <- ggplot(loadings, aes(Site, PC4)) +
   geom_point(aes(colour = factor(Chromosome))) +
   ggtitle("Plot of PC4 by genomic position in Z014 (B73xKy21)") +
   xlab("Genomic Position")

# PC5
p5 <- ggplot(loadings, aes(Site, PC5)) +
   geom_point(aes(colour = factor(Chromosome))) +
   ggtitle("Plot of PC5 by genomic position in Z014 (B73xKy21)") +
   xlab("Genomic Position")

library(ggpubr)
ggarrange(p1,p2, p3,p4,p5 + rremove("x.text"),
          ncol = 1, nrow = 5)
ggsave("first5PCs.png", plot = last_plot(), device = NULL, path = NULL)


# ----------------------------------------
# Try different genotype table
# http://datacommons.cyverse.org/browse/iplant/home/shared/panzea/genotypes/SNPs/NAM_map_and_genos-120731.zip
# 1106 markers from original NAM paper
# Doing PCA on markers for one nam family
# ky21 --> Z014
# ----------------------------------------

nam_map <- read.csv("ZeaGBSv27_publicSamples_imputedV5_AGPv4-181023_NAM_Filter.csv", header = FALSE)

# Only get info for one nam family (Z014)
nam_map$family <- nam_map[,1]
nam_map$family <- gsub("E[0-9]{4}", "", nam_map$family)
ky21_genos <- nam_map[which(nam_map$family == "Z014"), ]

# Add back in chromosome information
ky21_genos <- rbind(nam_map[1:3,], ky21_genos)

# transpose data
ky21_genos <- as.data.frame(t(ky21_genos))

# Change column names
colnames(ky21_genos) <- as.character(unlist((ky21_genos[1,])))
ky21_genos <- ky21_genos[-1,]

# Get values only
ky21_geno_only <- ky21_genos[-1107,4:length(ky21_genos)]

# Turn genotypes into values
ky21_geno_only <- as.data.frame(lapply(ky21_geno_only, function(x) as.numeric(as.character(x))))

# Run pca
colnames(ky21_geno_only) <- NULL
pca_ky21 <- prcomp(t(ky21_geno_only), scale = FALSE)

# Extract out PCs
loadings_1106Markers <- as.data.frame(pca_ky21$rotation)

# Add in site information
loadings_1106Markers <- cbind(ky21_genos[-1107,1:3],loadings_1106Markers)

# Add in a site to plot by
loadings_1106Markers$Site <- seq(1:nrow(loadings_1106Markers))
loadings_1106Markers$Chr <- as.numeric(as.character(loadings_1106Markers$Chr))

var <- cumsum(pca_ky21$sdev^2/sum(pca_ky21$sdev^2))

# PC1
p1 <- ggplot(loadings_1106Markers, aes(Site, PC1)) +
  geom_point(aes(colour = factor(Chr))) +
  ggtitle("Plot of PC1 by genomic position in Z014 (B73xKy21)") +
  xlab("Genomic Position") +
   ylim(-0.1, 0.1) +
   theme(axis.text.x=element_blank())

# PC2
p2 <- ggplot(loadings_1106Markers, aes(Site, PC2)) +
  geom_point(aes(colour = factor(Chr))) +
  ggtitle("Plot of PC2 by genomic position in Z014 (B73xKy21)") +
  xlab("Genomic Position")+
   ylim(-0.1, 0.1)+
   theme(axis.text.x=element_blank())

# PC3
p3 <- ggplot(loadings_1106Markers, aes(Site, PC3)) +
  geom_point(aes(colour = factor(Chr))) +
  ggtitle("Plot of PC3 by genomic position in Z014 (B73xKy21)") +
  xlab("Genomic Position")+
   ylim(-0.1, 0.1)+
   theme(axis.text.x=element_blank())

# PC4
p4 <- ggplot(loadings_1106Markers, aes(Site, PC4)) +
  geom_point(aes(colour = factor(Chr))) +
  ggtitle("Plot of PC4 by genomic position in Z014 (B73xKy21)") +
  xlab("Genomic Position")+
   ylim(-0.1, 0.1)+
   theme(axis.text.x=element_blank())

# PC5
p5 <- ggplot(loadings_1106Markers, aes(Site, PC5)) +
  geom_point(aes(colour = factor(Chr))) +
  ggtitle("Plot of PC5 by genomic position in Z014 (B73xKy21)") +
  xlab("Genomic Position")+
   ylim(-0.1, 0.1)+
   theme(axis.text.x=element_blank())

library(ggpubr)
ggarrange(p1,p2, p3,p4,p5 + rremove("x.text"),
          ncol = 1, nrow = 5)
ggsave("first5PCs_1106Markers.png", plot = last_plot(), device = NULL, path = NULL,
       width = 18, height = 11)


# ----------------------------------------
# PCA on genotypes/markers
# http://datacommons.cyverse.org/browse/iplant/home/shared/panzea/genotypes/SNPs/NAM_map_and_genos-120731.zip
# 1106 markers from original NAM paper
# Doing PCA on markers for all NAM families
# ----------------------------------------

# All nam familes, 1106 markers original NAM markers
nam_map <- read.csv("ZeaGBSv27_publicSamples_imputedV5_AGPv4-181023_NAM_Filter.csv", header = FALSE)

# Don't really need this anymore, keeping only becuase
# I don't want to change old code
# What it does: Only get info for one nam family (Z014)
nam_map$family <- nam_map[,1]
nam_map$family <- gsub("E[0-9]{4}", "", nam_map$family)

# transpose data
all_namGenos <- as.data.frame(t(nam_map))

# Change column names
colnames(all_namGenos) <- as.character(unlist((all_namGenos[1,])))
all_namGenos <- all_namGenos[-1,]

# Get values only
nam_geno_only <- all_namGenos[-1107,4:length(all_namGenos)]

# Turn genotypes into values
nam_geno_only <- as.data.frame(lapply(nam_geno_only, function(x) as.numeric(as.character(x))))

# Run pca
colnames(nam_geno_only) <- NULL
pca_nam <- prcomp(t(nam_geno_only), scale = FALSE)

# Extract out PCs
allNAM_loadings_1106Markers <- as.data.frame(pca_nam$rotation)

# Add in site information
allNAM_loadings_1106Markers <- cbind(all_namGenos[-1107,1:3], allNAM_loadings_1106Markers)

# Add in a site to plot by
allNAM_loadings_1106Markers$Site <- seq(1:nrow(allNAM_loadings_1106Markers))
allNAM_loadings_1106Markers$Chr <- as.numeric(as.character(allNAM_loadings_1106Markers$Chr))

# Cumulative sum of all PCs
var <- cumsum(pca_nam$sdev^2/sum(pca_nam$sdev^2))
plot(var, ylab = "Cumulative variance explained",
     xlab = "PCs")

# Plotting
# PC1
p1 <- ggplot(allNAM_loadings_1106Markers, aes(Site, PC1)) +
  geom_point(aes(colour = factor(Chr))) +
  ggtitle("Plot of PC1 by genomic position in all NAM") +
  xlab("Genomic Position") +
  ylim(-0.15, 0.15) +
  theme(axis.text.x=element_blank())

# PC2
p2 <- ggplot(allNAM_loadings_1106Markers, aes(Site, PC2)) +
  geom_point(aes(colour = factor(Chr))) +
  ggtitle("Plot of PC2 by genomic position in all NAM") +
  xlab("Genomic Position")+
  ylim(-0.15, 0.15)+
  theme(axis.text.x=element_blank())

# PC3
p3 <- ggplot(allNAM_loadings_1106Markers, aes(Site, PC3)) +
  geom_point(aes(colour = factor(Chr))) +
  ggtitle("Plot of PC3 by genomic position in all NAM") +
  xlab("Genomic Position")+
  ylim(-0.15, 0.15)+
  theme(axis.text.x=element_blank())

# PC4
p4 <- ggplot(allNAM_loadings_1106Markers, aes(Site, PC4)) +
  geom_point(aes(colour = factor(Chr))) +
  ggtitle("Plot of PC4 by genomic position in all NAM") +
  xlab("Genomic Position")+
  ylim(-0.15, 0.15)+
  theme(axis.text.x=element_blank())

# PC5
p5 <- ggplot(allNAM_loadings_1106Markers, aes(Site, PC5)) +
  geom_point(aes(colour = factor(Chr))) +
  ggtitle("Plot of PC5 by genomic position in all NAM") +
  xlab("Genomic Position")+
  ylim(-0.15, 0.15)+
  theme(axis.text.x=element_blank())

# PC6
p6 <- ggplot(allNAM_loadings_1106Markers, aes(Site, PC6)) +
  geom_point(aes(colour = factor(Chr))) +
  ggtitle("Plot of PC6 by genomic position in all NAM") +
  xlab("Genomic Position")+
  ylim(-0.15, 0.15)+
  theme(axis.text.x=element_blank())

# PC7
p7 <- ggplot(allNAM_loadings_1106Markers, aes(Site, PC7)) +
  geom_point(aes(colour = factor(Chr))) +
  ggtitle("Plot of PC7 by genomic position in all NAM") +
  xlab("Genomic Position")+
  ylim(-0.15, 0.15)+
  theme(axis.text.x=element_blank())

# PC8
p8 <- ggplot(allNAM_loadings_1106Markers, aes(Site, PC8)) +
  geom_point(aes(colour = factor(Chr))) +
  ggtitle("Plot of PC8 by genomic position in all NAM") +
  xlab("Genomic Position")+
  ylim(-0.15, 0.15)+
  theme(axis.text.x=element_blank())

# PC9
p9 <- ggplot(allNAM_loadings_1106Markers, aes(Site, PC9)) +
  geom_point(aes(colour = factor(Chr))) +
  ggtitle("Plot of PC9 by genomic position in all NAM") +
  xlab("Genomic Position")+
  ylim(-0.15, 0.15)+
  theme(axis.text.x=element_blank())

# PC10
p10 <- ggplot(allNAM_loadings_1106Markers, aes(Site, PC10)) +
  geom_point(aes(colour = factor(Chr))) +
  ggtitle("Plot of PC10 by genomic position in all NAM") +
  xlab("Genomic Position")+
  ylim(-0.15, 0.15)+
  theme(axis.text.x=element_blank())

library(ggpubr)
ggarrange(p1,p2, p3,p4,p5,p6,p7,p8,p9,p10 + rremove("x.text"),
          ncol = 1, nrow = 10)
ggsave("NAM_wide_first10PCs_1106Markers.png", plot = last_plot(), device = NULL, path = NULL,
       width = 15, height = 22)


# ----------------------------------------
# Try different genotype table
# http://datacommons.cyverse.org/browse/iplant/home/shared/panzea/genotypes/GBS/v23/NAM.GBSv23.imputedMarkers.0.2cm.zip
# /iplant/home/shared/panzea/genotypes/GBS/v23/NAM.GBSv23.imputedMarkers.0.2cm.zip
#  markers from Eli Rodgers-Melnick and 
# Peter Bradbury imputation
# Doing PCA on markers for all NAM families
# ----------------------------------------

# Change directory
setwd("/Users/mbb262-admin/Box\ Sync/Cornell_PhD/labProjects/hap_gwa/pcs_fastAssociation_282/NAM.GBSv23.imputedMarkers.0.2cm")

# All files
files <- list.files()

# Read in data
combined_files <- do.call("rbind", lapply(files, read.table, header = TRUE))

# Get values only
allNAM_7kmarkers <- combined_files[,6:ncol(combined_files)]

# Run pca
colnames(allNAM_7kmarkers) <- NULL
pca_nam7k <- prcomp(t(allNAM_7kmarkers), scale = FALSE)

# Extract out PCs
# Takes a long time
allNAM_loadings_7kmarkers <- as.data.frame(pca_nam7k$rotation)

# Add in site information
allNAM_loadings_7kmarkers <- cbind(combined_files[,3:4], allNAM_loadings_7kmarkers)

# Sort chromosomes first (10 still coming after 1)
allNAM_loadings_7kmarkers <- allNAM_loadings_7kmarkers[order(allNAM_loadings_7kmarkers$chr),] 


# Add in a site to plot by
allNAM_loadings_7kmarkers$Site <- seq(1:nrow(allNAM_loadings_7kmarkers))
allNAM_loadings_7kmarkers$chr <- as.numeric(as.character(allNAM_loadings_7kmarkers$chr))

# Cumulative sum of all PCs
var <- cumsum(pca_nam7k$sdev^2/sum(pca_nam7k$sdev^2))
plot(var, ylab = "Cumulative variance explained",
     xlab = "PCs")
tmp <- data.frame(var)

# Plotting
library(ggplot2)

# PC1
p1 <- ggplot(allNAM_loadings_7kmarkers, aes(Site, PC1)) +
  geom_point(aes(colour = factor(chr))) +
  ggtitle("Plot of PC1 by genomic position in all NAM (7k markers)") +
  xlab("Genomic Position") +
  ylim(-0.05, 0.05) +
  theme(axis.text.x=element_blank())

# PC2
p2 <- ggplot(allNAM_loadings_7kmarkers, aes(Site, PC2)) +
  geom_point(aes(colour = factor(chr))) +
  ggtitle("Plot of PC2 by genomic position in all NAM (7k markers)") +
  xlab("Genomic Position")+
  ylim(-0.05, 0.05)+
  theme(axis.text.x=element_blank())

# PC3
p3 <- ggplot(allNAM_loadings_7kmarkers, aes(Site, PC3)) +
  geom_point(aes(colour = factor(chr))) +
  ggtitle("Plot of PC3 by genomic position in all NAM (7k markers)") +
  xlab("Genomic Position")+
  ylim(-0.05, 0.05)+
  theme(axis.text.x=element_blank())

# PC4
p4 <- ggplot(allNAM_loadings_7kmarkers, aes(Site, PC4)) +
  geom_point(aes(colour = factor(chr))) +
  ggtitle("Plot of PC4 by genomic position in all NAM (7k markers)") +
  xlab("Genomic Position")+
  ylim(-0.05, 0.05)+
  theme(axis.text.x=element_blank())

# PC5
p5 <- ggplot(allNAM_loadings_7kmarkers, aes(Site, PC5)) +
  geom_point(aes(colour = factor(chr))) +
  ggtitle("Plot of PC5 by genomic position in all NAM (7k markers)") +
  xlab("Genomic Position")+
  ylim(-0.05, 0.05)+
  theme(axis.text.x=element_blank())

# PC6
p6 <- ggplot(allNAM_loadings_7kmarkers, aes(Site, PC6)) +
  geom_point(aes(colour = factor(chr))) +
  ggtitle("Plot of PC6 by genomic position in all NAM (7k markers)") +
  xlab("Genomic Position")+
  ylim(-0.05, 0.05)+
  theme(axis.text.x=element_blank())

# PC7
p7 <- ggplot(allNAM_loadings_7kmarkers, aes(Site, PC7)) +
  geom_point(aes(colour = factor(chr))) +
  ggtitle("Plot of PC7 by genomic position in all NAM (7k markers)") +
  xlab("Genomic Position")+
  ylim(-0.05, 0.05)+
  theme(axis.text.x=element_blank())

# PC8
p8 <- ggplot(allNAM_loadings_7kmarkers, aes(Site, PC8)) +
  geom_point(aes(colour = factor(chr))) +
  ggtitle("Plot of PC8 by genomic position in all NAM (7k markers)") +
  xlab("Genomic Position")+
  ylim(-0.05, 0.05)+
  theme(axis.text.x=element_blank())

# PC9
p9 <- ggplot(allNAM_loadings_7kmarkers, aes(Site, PC9)) +
  geom_point(aes(colour = factor(chr))) +
  ggtitle("Plot of PC9 by genomic position in all NAM (7k markers)") +
  xlab("Genomic Position")+
  ylim(-0.05, 0.05)+
  theme(axis.text.x=element_blank())

# PC10
p10 <- ggplot(allNAM_loadings_7kmarkers, aes(Site, PC10)) +
  geom_point(aes(colour = factor(chr))) +
  ggtitle("Plot of PC10 by genomic position in all NAM (7k markers)") +
  xlab("Genomic Position")+
  ylim(-0.05, 0.05)+
  theme(axis.text.x=element_blank())

library(ggpubr)
ggarrange(p1,p2, p3,p4,p5,p6,p7,p8,p9,p10 + rremove("x.text"),
          ncol = 1, nrow = 10)
ggsave("NAM_wide_first10PCs_7k_Markers.png", plot = last_plot(), device = NULL, path = NULL,
       width = 16, height = 22)


# -----------------------------------------
# Plotting one nam family by top 5 PCs
# PCA on Genotypes
# `/iplant/home/shared/panzea/genotypes/GBS/v27/ZeaGBSv27_publicSamples_imputedV5_AGPv4-181023.vcf.gz`
# First dataset again for only ky21
# Sataset ran through TASSEL function Numerical Genotype
# aka collapse
# ----------------------------------------

meanImpute <- function(X) {
  apply(X, 2, function(x) {
    x[is.na(x)] <- mean(x, na.rm=TRUE)
    return(x)
  })
}

# Change directory
setwd("/Users/mbb262-admin/Box\ Sync/Cornell_PhD/labProjects/hap_gwa/pcs_fastAssociation_282")

# read genotype table
geno <- read.table("ZeaGBSv27_publicSamples_imputedV5_AGPv4-181023_NAM_Filter_MAF0.35_with_Probability.txt", header = T)

# Genotyoes (only numeric)
geno2 <- geno[,2:length(geno)]

# Turn genotypes into values (recheck this, think I'm setting values strangely)
# geno3 <- as.data.frame(lapply(geno2, function(x) as.numeric(x)))

# Transpose markers
geno3 <- t(geno2)

# Run PCA on markers
results_pca <- prcomp(t(na.omit(geno3)), scale = FALSE)

# Extract out PCs
loadings <- as.data.frame(results_pca$rotation)

# Add in site information
# siteNames <- read.table("position_list_ky21.txt", header = TRUE)
tmp <- t(na.omit(geno))
loadings <- cbind(tmp,loadings)

# Plot pcas by site
plot(loadings$tmp, abs(loadings$PC1))

# ggplot, color by chromosome
library(ggplot2)

# PC1
p1 <- ggplot(loadings, aes(Site, abs(PC1))) +
  geom_point(aes(colour = factor(Chromosome))) +
  ggtitle("Plot of PC1 by genomic position in Z014 (B73xKy21)") +
  xlab("Genomic Position")



