# ---------------------------------------
# Merritt Burch
# mbb262@cornell.edu
# 2019-08-26
# Script to look at a distance matrix 
# produced by TASSEL and filter on relatedness
# ---------------------------------------

# Set workdir
setwd("~/git_projects/haplotype_v_snp_gwas/src/main/shell")

# Read in data
distance <- read.table("ZeaGBSv27_20171204_AGPv4_Ames_genotype_distance4.txt", header = FALSE)

# Turn matrix values into numeric
numeric_distance <- as.data.frame(lapply(distance[,2:ncol(distance)], function(x) as.numeric(as.character(x))))

# Turn into matrix
temp <- as.matrix(numeric_distance)

# Plot whole matrix
hist(temp)
hist(temp, breaks = 200, xlim = c(0,0.4))


# -----------------------------
# Subset matrix on relatedness 
# -----------------------------

# Hide upper triangle
upper <- distance
upper[upper.tri(distance)] <- NA
upper<-as.data.frame(upper)

# Plot whole matrix
hist(as.matrix(upper[,2:ncol(upper)]),
     xlab = "Distance matrix values in Ames (1 - IBS (identity by state) similarity)",
     main = "Distance matrix similarity for all Ames individuals and markers")

hist(as.matrix(upper[,2:ncol(upper)]), breaks = 200, xlim = c(0,0.4),
     xlab = "Distance matrix values in Ames (1 - IBS (identity by state) similarity)",
     main = "Distance matrix similarity for all Ames individuals and markers (Zoomed in)")

# Filter bottom part of distance matrix for individuals
# related less than 25% (may adjust later)
# Gives the lines with >0.25 relatedness, take list and subset main distaance file
library(dplyr)
ames_lines_subset_relatedness <- upper %>%
  filter_all(any_vars(. > 0.25)) 

# Extract taxa IDs and filter from main distance file
highly_related_ames_taxa <- data.frame(ames_lines_subset_relatedness[,1])
all_ames_taxa <- distance
# as.data.frame(distance$V1)

# Subset all taxa from distance matrix by the highly related taxa
# aka: remove highly related taaxa IDs from all taxa IDs to be left
# with the not-very-related taxa
# Should be expecting a df with 3545-251 = 3294 rows (it works!)
not_very_related_ames  <- data.frame(all_ames_taxa[!(all_ames_taxa[,1] %in% highly_related_ames_taxa[,1]),])

# Check with a histogram
# It looks like it worked
hist(as.matrix(not_very_related_ames[,2:ncol(not_very_related_ames)]), breaks = 200)

# Clean it up
temp <- data.frame(not_very_related_ames[,1])
colnames(temp) <- "Taxa"

# Save taxa list
write.table(temp, file = "Ames_taxa_lessThan0.25_related.txt",
            row.names = FALSE, quote = FALSE)


# ------------------
# K-means clustering
# ------------------

# Try Anju's method

# Load distance matrix
distance <- read.table("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/ames_tests/ZeaGBSv27_20171204_AGPv4_Ames_genotype_distance4.txt", header = FALSE)

# Distance matrix on windows machine
distance <- read.table("C:/Users/merri/Box Sync/Cornell_PhD/labProjects/hap_gwa/ames_tests/ZeaGBSv27_20171204_AGPv4_Ames_genotype_distance4.txt", header = FALSE)

# New distance matrix
distance <- read.table("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/ames_tests/MatrixAmes_MinCount500_MAF0.05_removedMissing.txt", header = FALSE)

#packages to load
library(factoextra)
library(cluster)
library(NbClust)

#number of centroids
k=150 

# Subset data and omit missing data
rownames(distance) <- distance$V1
dat <- na.omit(distance[,2:ncol(distance)])

# Run clustering
model <- kmeans(dat,k)

# Place to store all of the taxa
index <- c()

#calculate indices of closest instance to centroid
for (i in 1:k){
  rowsum <- rowSums(abs(dat[which(model$cluster==i),] - model$centers[i,]))
  index[i] <- names(which.min(rowsum))
}

# Turn into dataframe and fix colnames
index <- as.data.frame(index)
colnames(index) <- "Taxa"

# Save taxa in txt file
write.table(index, "k_means_150clusters_Ames_taxa.txt", row.names = FALSE, quote = FALSE)


# ---------------------
# Do PCA on these taxa
# ---------------------

# Set directory
setwd("C:/Users/merri/Box Sync/Cornell_PhD/labProjects/hap_gwa/ames_tests")

# read genotype table
geno <- read.table("ZeaGBSv27_20171204_AGPv4_Ames_Kmeans150_MinCount60_MAF0.05_KNNimp_with_Probability_numeric.txt", header = T)

# Genotyoes (only numeric)
geno2 <- geno[,2:length(geno)]

# Transpose markers
geno3 <- t(geno2)

# Run PCA on markers
results_pca <- prcomp(t(na.omit(geno3)), scale = FALSE)

# Extract out PCs
loadings <- as.data.frame(results_pca$rotation)

# Add in site information
# siteNames <- read.table("position_list_ky21.txt", header = TRUE)
tmp <- t(na.omit(geno3))
loadings <- cbind(tmp,loadings)

# Plot pcas by site
plot(loadings$tmp, abs(loadings$PC1))


# Get relevant Ames phenos
ames_phenos <- read.table("C:/Users/merri/git_projects/haplotype_v_snp_gwas/data/all_Ames_Phenos.txt",header = TRUE)
ames_subset <- ames_phenos[,c(1,6:14)]

# Add in long names (GBS format)
long_names <- data.frame(distance$V1)
long_names$copy <- distance$V1
long_names$copy <- gsub(":[0-9]{9}", "", long_names$copy)
long_names$copy <- tolower(long_names$copy)

# Merge datasets by names
together <- merge(x = ames_subset, y = long_names, by.x = "Entry_ID", by.y = "copy")

# Change order
together <- together[,c(11,2:10)]

# Save
write.table(together, "ames_phenos_subset.txt", row.names = FALSE, quote = FALSE, sep = "\t")


# ------------------
# Plot gWAS results from
# using 60, 10, and 5 PCs
# to control for populations str
# -------------------------
setwd("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/ames_tests")

ames60 <- read.table("FastAssociation_PC_Ames_Kmeans150_MinCount60_MAF0.05_KNNimp_60PCs + Filtered_BLUEs_ames_phenos_GDD_DTS + Ames_MinCount500_MAF0.05_removedMissing.txt",
                     header = TRUE)

ames10 <- read.table("FastAssociation_PC_ZeaGBSv27_20171204_AGPv4_Ames_Kmeans150_MinCount60_MAF0.05_KNNimp_10PCs + Filtered_BLUEs_ames_phenos_GDD_DTS + Ames_MinCount500_MAF0.05_removedMissing.txt",
                     header = TRUE)

ames5 <- read.table("FastAssociation_PC_Ames_Kmeans150_MinCount60_MAF0.05_KNNimp_5PCs + Filtered_BLUEs_ames_phenos_GDD_DTS + Ames_MinCount500_MAF0.05_removedMissing.txt",
                    header = TRUE)


# Plot results

# ggman
library(ggman)
library(qqman)

# 60
png("GDD_DTA_Ames_60.png",width = 5, height = 3.5, units = "in", res = 500)
ggman(ames60, chrom="Chr", bp="Pos", snp="Marker", pvalue="p", pointSize = 0.3, lineColour = "black",
      title = "GDD_DTS in Ames (Romay 2013) + 60 PCs", sigLine = -log10(0.05))
dev.off()
#10
png("GDD_DTA_Ames_10.png",width = 5, height = 3.5, units = "in", res = 500)
ggman(ames10, chrom="Chr", bp="Pos", snp="Marker", pvalue="p", pointSize = 0.3, lineColour = "black",
      title = "GDD_DTS in Ames (Romay 2013) + 10 PCs", sigLine = -log10(0.05))
dev.off()
#5
png("GDD_DTA_Ames_5.png",width = 5, height = 3.5, units = "in", res = 500)
ggman(ames5, chrom="Chr", bp="Pos", snp="Marker", pvalue="p", pointSize = 0.3, lineColour = "black",
      title = "GDD_DTS in Ames (Romay 2013) + 5 PCs", sigLine = -log10(0.05))
dev.off()


# Anovas on flowering time and PCs

pcs_phenos <- read.table("PC_Ames_Kmeans150_MinCount60_MAF0.05_KNNimp_60PCsv_2.txt",
                         header = TRUE)
tenPCs <- pcs_phenos[,c(1:11,62)]
model_5PCs_DTS <- lm(GDD_DTS_BLUP_Romay2013~PC1+PC2+PC3+PC4+PC5, tenPCs)
summary(model_5PCs_DTS)
plot(model_5PCs_DTS)


model_10PCs_DTS <- lm(GDD_DTS_BLUP_Romay2013~PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, tenPCs)
summary(model_10PCs_DTS)
plot(model_10PCs_DTS)







