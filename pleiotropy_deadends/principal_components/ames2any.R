# --------------------------------------------
# Merritt Burch
# mbb262@cornell.edu
# 2019-11-06
#
# Script to do local PCA in Ames and apply
# loadings/rotations to any other pop
# --------------------------------------------

# Load packages
library(gdsfmt)
library(SNPRelate)

# Load whole Beagle imputed SNPs
vcf.fn <- "~/Box Sync/Cornell_PhD/labProjects/hackathons/2019-10-21_hackathon/ZeaGBSv27_20171204_AGPv4_Ames_filteredGenotypes_beagleImputation_Filter_intersectSitesWithNAM_Filter.vcf"
beagle_vcf <- snpgdsVCF2GDS(vcf.fn, "beagle_filtered.gds", method="copy.num.of.ref", snpfirstdim=TRUE, ignore.chr.prefix = "S")
genofile_ames <- snpgdsOpen(beagle_vcf)

# Load in subset of snps
vcf.fn <- "~/Box Sync/Cornell_PhD/labProjects/hackathons/2019-10-21_hackathon/ZeaGBSv27_20171204_AGPv4_Ames_filteredGenotypes_beagleImputation_Filter_intersectSitesWithNAM_Filter_Filter-smallersubset.vcf"
beagle_vcf <- snpgdsVCF2GDS(vcf.fn, "beagle_filtered.gds", method="copy.num.of.ref", snpfirstdim=TRUE, ignore.chr.prefix = "S")
genofile_ames <- snpgdsOpen(beagle_vcf)


# -----------------
# Useful functions
# -----------------

# Function to get Ames SNP means and V matrix from global PCs
get_ames_means_Vmatrix <- function(genofileID, num_gPCs){
  
  # Turn genotypes into numeric values
  X_list <- snpgdsGetGeno(genofileID, snpfirstdim=FALSE, with.id=TRUE, verbose = FALSE)
  
  # Get SNP means
  centers <- colMeans(X_list$genotype)
  
  # Do PCA
  V_global <- prcomp(X_list$genotype, center = TRUE, rank. = num_gPCs)$rotation
  
  # Return V_global matrix and Ames SNP means
  return(list(centers, V_global))
}


# Function to apply Ames rotations to any other population
ames2any <- function(genofileID, snpMeans, Vglobal){
  
  # Turn genotypes into numeric values
  X_list <- snpgdsGetGeno(genofileID, snpfirstdim=FALSE, with.id=TRUE, verbose = FALSE)
  
  # Scale matrix
  PC_adj <- scale(X_list$genotype, center = snpMeans, scale = FALSE) %*% Vglobal
  
  # Return scaled matrix
  return(PC_adj)
}


# Test the function on Ames data
ames_means_V <- get_ames_means_Vmatrix(genofile_ames, 3) # Takes 20 mins to run

# Adjust and create PCs
pc_adj <- ames2any(genofile_ames, snpMeans = ames_means_V[[1]], Vglobal = ames_means_V[[2]])


# ----------------------------------------
# Apply function to a different population
# ----------------------------------------

# Load  data
nam.vcf <- "~/Box Sync/Cornell_PhD/labProjects/hackathons/2019-10-21_hackathon/intersectJoin_NAMbyAmes_shared_beagleImputed.vcf.gz"
nam_vcf <- snpgdsVCF2GDS(nam.vcf, "nam.gds", method="copy.num.of.ref", snpfirstdim=TRUE, ignore.chr.prefix = "S")
genofile_nam <- snpgdsOpen(nam_vcf)

# Use function
pc_adj_nam <- ames2any(genofile_nam, snpMeans = ames_means_V[[1]], Vglobal = ames_means_V[[3]])

# Plot out PC axes
plot(ames_means_V[[2]]$x[,1], ames_means_V[[2]]$x[,2])
plot(pc_adj_nam[,1], pc_adj_nam[,2])

# Plot with colors
nam_subpops <- read.csv("~/git_projects/haplotype_v_snp_gwas/data/nam_subpops.csv")
ames_nam_color_sub <- data.frame(cbind(read.gdsn(index.gdsn(genofile_nam, "sample.id")), pc_adj_nam))
ames_nam_color_sub[,1] <- gsub("E[0-9]{4}:[0-9]{9}", "", ames_nam_color_sub[,1])
ames_nam_color_sub$PC1 <- as.numeric(as.character(ames_nam_color_sub$PC1))
ames_nam_color_sub$PC2 <- as.numeric(as.character(ames_nam_color_sub$PC2))
ames_nam_color_sub <- merge(x = ames_nam_color_sub, y = nam_subpops, by.x = "V1", by.y = "taxa")
library(ggplot2)
ggplot(ames_nam_color_sub, aes(x=PC1, y=PC2, color = V1)) +
  geom_point(aes(shape = subpop))


# NAM only PCs 
nam_numeric <- snpgdsGetGeno(genofile_nam, snpfirstdim=FALSE, with.id=TRUE, verbose = FALSE)
pca_nam <- prcomp(nam_numeric$genotype, center = TRUE, rank. = 3)
plot(pca_nam$x[,1], pca_nam$x[,2])

# Color stuff
nam_color <- data.frame(cbind(read.gdsn(index.gdsn(genofile_nam, "sample.id")), pca_nam$x))
nam_color_sub <- merge(x = nam_color, y = nam_subpops, by.x = "V1", by.y = "taxa")
nam_color[,1] <- gsub("E[0-9]{4}:[0-9]{9}", "", nam_color[,1])
nam_color$PC1 <- as.numeric(as.character(nam_color$PC1))
nam_color$PC2 <- as.numeric(as.character(nam_color$PC2))
ggplot(nam_color_sub, aes(x=PC1, y=PC2, color = V1)) +
  geom_point(aes(shape = subpop))