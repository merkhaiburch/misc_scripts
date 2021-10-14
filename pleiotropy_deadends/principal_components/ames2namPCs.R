# ---------------------------------------
# Merritt Burch
# mbb262@cornell.edu
# 2019-10-28
# Script to apply the coefficients of the
# local PCs from the Ames populations to 
# the NAM population
#
# Then extract local NAM PCs and apply them in
# the models:

# y ~ globalPCs+localPCs (PCs from NAM only)
# y ~ SNPs(400k) + globalPCs + localPCs (PCs from NAM only)

# y ~ (globalPCs + localPCs from ames applied to 40k NAM SNPs)

# Where y = DTS from Buckler 2009
# ---------------------------------------

# Steps to running local PCA on a genotype (X) matrix
# - Import vcf file
# - get genotype (X) matrix (SNPRelate gives count of reference allele, normally programs giive alterate counts)
# - Scale X matrix with varable number of PCs (1-5)
# - Calculate local PC windows (using list of SNP IDs) on scaled matrix
# - Extract out 1-3 local PC loadings from each analysis, add to one central list
# - Apply loadings from Ames SNPs to NAM SNPs 

# -------------------------------------------
# Calculate/adjust local PCs with global PCs
# -------------------------------------------

# Load packages
library(SNPRelate)
library(gdsfmt)
library(ggplot2)

# Set directory
setwd("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/ames_tests/Sep25_localPCs")

# Ames and NAM sites should be exactly the same between these two samples

# Load in Beagle imputed SNPs that were filtered with the 'temp' (aka ames_taxa_gbs_deduplicated.txt, code@bottom)
vcf.fn <- "~/Box Sync/Cornell_PhD/labProjects/hackathons/2019-10-21_hackathon/ZeaGBSv27_20171204_AGPv4_Ames_filteredGenotypes_beagleImputation_Filter_intersectSitesWithNAM_Filter.vcf"
beagle_vcf <- snpgdsVCF2GDS(vcf.fn, "beagle_filtered.gds", method="copy.num.of.ref", snpfirstdim=TRUE, ignore.chr.prefix = "S")
snpgdsSummary("beagle_filtered.gds")
genofile_ames <- snpgdsOpen(beagle_vcf)

# Load in NAM SNPs that are shared (intersect) with the Ames SNPs (did intersection in TASSEL)
nam.vcf <- "~/Box Sync/Cornell_PhD/labProjects/hackathons/2019-10-21_hackathon/intersectJoin_NAMbyAmes_shared_beagleImputed.vcf.gz"
nam_vcf <- snpgdsVCF2GDS(nam.vcf, "nam.gds", method="copy.num.of.ref", snpfirstdim=TRUE, ignore.chr.prefix = "S")
snpgdsSummary("nam.gds")
genofile_nam <- snpgdsOpen(nam_vcf)


# -------------------------------
# Calculate global and local PCs
# -------------------------------

# Run global function from local_pca_ames.R
# Centers and decomposes matrix
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

# Calculate PCs this way --> centered only om Ames
G <- snpgdsGRM(genofile_ames, method="EIGMIX") ## centered only
pca_global <- decomp(G$grm, k=5, center = TRUE)


# Local PC function the same as local_pca_ames.R EXCEPT 
# that it returns loadings (aka rotations) in matrices within a list

# Function to adjust a genotype matrix with covariates
regress_out <- function(X, Q, include_intercept=TRUE) {
  
  n <- nrow(X)
  
  stopifnot(n == nrow(Q))

  # Include intercept to Q
  if (include_intercept & sum(apply(X, 2, var) ==0) == 0) {
    Q <- cbind(1, Q)
  }
  
  # Matrix of projection
  H <- diag(n) - Q %*% solve(crossprod(Q)) %*% t(Q)
  
  # Projecting X out of Q
  H %*% X
}

# Collector of numbers
all_local_pcs <- list()

# The function to calculate local PCs
calc_local_pcs_returnRotations <- function(genofileID, globalPCs, numLocalPcsOut, windowSize = 1e7){
  # Iterate over chromosmes
  for (chrom in seq(1:10)){
    
    # Give a status message
    # cat("\n I am analyzing chromosome", chrom, "\n")
    
    # Get whole SNP list
    snp_map <- snpgdsSNPList(genofileID)
    
    # Subset SNPs by chromosome and window
    windows_genome <- cut_width(snp_map$position[snp_map$chromosome == as.character(chrom)], width = windowSize)
    
    # Subset SNPs by chromosome, window, put them into separate lists to iterate over
    snp_split <- split(x=snp_map$snp.id[snp_map$chromosome == as.character(chrom)], 
                       f=windows_genome)
    
    # Another loop to iterate over all windows and calculate local PCs
    collector <- c()
    for (windows in seq(1:length(snp_split))){
      #length(snp_split)
      # Print status
      # cat("\n I am on window ", windows, "/", length(snp_split), "\n")
       
      # Get genotype matrix
      # Genotype data --> get X matrix
      ## snpfirstdim=FALSE to get a nxp matrix
      ## with.id=TRUE to get sample and SNP information in a list
      X_list <- snpgdsGetGeno(genofileID, snpfirstdim=FALSE, with.id=TRUE, 
                              snp.id = unlist(snp_split[windows]), verbose = F)
      
      # Format X matrix
      X <- X_list$genotype
      rownames(X) <- X_list$sample.id
      colnames(X) <- X_list$snp.id
      
      ## to get alternate counts instead reference counts
      X <- 2-X
      
      # Adjust subset of snps using global PCs
      # No longer has taxa ids
      adj_snps <- regress_out(X, globalPCs)
      
      # Calculate and get local PCs
      # prcomp is ok with nxp matrix
      local_PCs <- prcomp(adj_snps)$rotation
      
      # Add local PCs to main dataframe
      all_local_pcs[[length(all_local_pcs)+1]] <- as.matrix(local_PCs[,1:numLocalPcsOut])
    }
  }
  # Return all the PCs in a big list
  return(all_local_pcs)
}

# Run the local PC function to get local PC loadings from ames
# Provides the [pxd] matrix
adj_pcs <- pca_global$PC[,1:5] # Global PCs to adjust by, 5 in this case
localPCs_ames2nam <- calc_local_pcs_returnRotations(genofile_ames, adj_pcs, 
                                                    numLocalPcsOut = 3, windowSize = 1e7)

# Example output
head(localPCs_ames2nam[[1]])


# --------------------------------
# Matrix multiply the whole thing
# --------------------------------

# For Ames we now have a [nxp] ("genofile_ames") and [pxd] ("localPCs_ames2nam", loadings) matrix
# For NAM we have a [nxp] matrix ("genofile_nam") with shared p between the two samples (the SNPs)

# Need to break up the NAM [nxp] matrix to be in consistent windows as Ames
# and matrix multiply the Ames loadings for the same SNPs

# Get list with matrices of NAM SNPs to adjust
# Gets local [nxp] matrices
all_nam_snps_windows <- list()
nam_snps_subsettedWindow <- function(genofileID, windowSize = 1e7){
  # Iterate over chromosmes
  for (chrom in seq(1:10)){
    # Get whole SNP list
    snp_map <- snpgdsSNPList(genofileID)
    
    # Subset SNPs by chromosome and window
    windows_genome <- cut_width(snp_map$position[snp_map$chromosome == as.character(chrom)], 
                                width = windowSize)
    
    # Subset SNPs by chromosome, window, put them into separate lists to iterate over
    snp_split <- split(x=snp_map$snp.id[snp_map$chromosome == as.character(chrom)], 
                       f=windows_genome)
    
    # Another loop to iterate over all windows and calculate local PCs
    collector <- c()
    
    for (windows in seq(1:length(snp_split))){
      
      # Get genotype (aka X) matrix
      ## snpfirstdim=FALSE to get a nxp matrix
      ## with.id=TRUE to get sample and SNP information in a list
      X_list <- snpgdsGetGeno(genofileID, snpfirstdim=FALSE, 
                              with.id=TRUE, 
                              snp.id = unlist(snp_split[windows]), 
                              verbose = F)
      
      # Format X matrix
      X <- X_list$genotype
      rownames(X) <- X_list$sample.id
      colnames(X) <- X_list$snp.id
      
      ## to get alternate counts instead reference counts
      X <- 2-X

      # Add local PCs to main dataframe
      all_nam_snps_windows[[length(all_nam_snps_windows)+1]] <- X
    }
  }
  # Return all the PCs in a big list
  return(all_nam_snps_windows)
}

# Run the local SNP window function
nam_snps_window_list <- nam_snps_subsettedWindow(genofile_nam, windowSize = 1e7)

# Create a function to iterate through both lists 
# to apply loadings from Ames SNPs to NAM SNPs within local windows
collect_adj_nam_snps <- list()

for (i in 1:length(localPCs_ames2nam)){
  collect_adj_nam_snps[[length(collect_adj_nam_snps)+1]] <- regress_out(t(nam_snps_window_list[[i]]), 
                                                                        localPCs_ames2nam[[i]], 
                                                                        include_intercept = F)
}

# Preview a few different [n2xd] matrices
data.frame(collect_adj_nam_snps[[2]])[1:10,1:5]

# Gather all SNPs into a single matrix
tmp_test <- lapply(collect_adj_nam_snps, t)
tmp_test_new <- do.call("cbind", tmp_test)

# Merge in phenotypic data
subset_nam_phenos <- all_NAM_phenos[,c(2,8:10)]
tmp_test_new <- cbind(row.names(tmp_test_new), data.frame(tmp_test_new))
tmp_test_new$copied_taxa <- gsub(":[0-9]{9}", "", tmp_test_new[,1])
namByAmes_snps_phenos <- merge(x = subset_nam_phenos, y = tmp_test_new,
                               by.x = "Geno_Code", by.y = "copied_taxa")

# Rearrange
namByAmes_snps_phenos_DTS <- namByAmes_snps_phenos[,c(1,3,6:ncol(namByAmes_snps_phenos))]

# Test SNPs against flowering time in model (GWAS)
dts_namByAmes <- lm(Days_To_Silk_BLUP_Sum0607_Buckler2009 ~ ., namByAmes_snps_phenos_DTS[,-1])


# ------------------------------------------------------
# Try calculating NAM specific local PCs
# Not using 'applied' NAM SNPs from the Ames population
# Calculate type III SS
# ------------------------------------------------------

# Import file with all PCs and GDD_DTS data
all_NAM_phenos <- read.csv("~/Box Sync/Cornell_PhD/labProjects/hackathons/2019-10-21_hackathon/all_NAM_phenos.txt", sep="")

# Calculate global PCs
G_nam <- snpgdsGRM(genofile_nam, method="EIGMIX") ## centered only
pca_global_nam <- decomp(G_nam$grm, k=5, center = TRUE)

# Calculate local PCs using the function from local_pca_ames.R
adj_pcs_nam <- pca_global_nam$PC[,1:5] # Global PCs to adjust by, 5 in this case
localPCs_1Mbp_nam <- calc_local_pcs(genofile_nam, adj_pcs_nam, numLocalPcsOut = 3, windowSize = 1e7)

# Combine global and local PCs, format taxa names
nam_global_local <- cbind(read.gdsn(index.gdsn(genofile_nam, "sample.id")), 
                          data.frame(pca_global_nam$PC[,1:5]), 
                          data.frame(localPCs_1Mbp_nam))
nam_global_local$copied_taxa <- gsub(":[0-9]{9}", "", nam_global_local[,1])

# Merge dataframes together
merged_pcs_phenos_nam <- merge(x = all_NAM_phenos, y = nam_global_local,
                           by.x = "Geno_Code", by.y = "copied_taxa")


# subset to only flowering time data and the PCs
subset_nam_FT_PCs <- merged_pcs_phenos_nam[,c(486,1,8:10,487:ncol(merged_pcs_phenos_nam))]

# Calculate type III sums of squares
# Set options/contrasts that sets sums to zero
options(contrasts = c("contr.sum", "contr.poly"))

# Store the model:
model_nam <- lm(Days_To_Silk_BLUP_Sum0607_Buckler2009 ~ ., subset_nam_FT_PCs[,-c(1:3,5)])

# Finally, call the drop1 function on each model component
# The results give the type III SS, including the p-values from an F-test.
# Not through car package, this is through base R
method1_baseR_typeIII <- drop1(model_nam, .~., test = "F")

# subset out 5 global PCs and intercept
method1_baseR_typeIII <- method1_baseR_typeIII[-c(1:6),]

# Make rownames one column
method1_baseR_typeIII <- data.table::setDT(method1_baseR_typeIII, keep.rownames = TRUE)[]

# Grep out chromosome, window, and PC info for plotting
method1_baseR_typeIII$chr <- gsub("_window_.*", "", method1_baseR_typeIII$rn)
method1_baseR_typeIII$window <- gsub(".*window_", "", method1_baseR_typeIII$rn)
method1_baseR_typeIII$PC <- gsub(".*_PC", "", method1_baseR_typeIII$window)
method1_baseR_typeIII$window <- as.numeric(gsub("_PC[0-9]", "", method1_baseR_typeIII$window))

# Plot in stacked barchart
ggplot(method1_baseR_typeIII, aes(x = window, y = `Sum of Sq`, fill = PC)) +
  geom_bar(stat = 'identity') +
  facet_grid(. ~ chr, scales = "free") +
  ggtitle("Type III Sums of Squares by top 3 local PCs from y(GDD_DTS) ~ 5 gPCs + 660 lPCs")


# ------------------------------------
# Get ready for tassel Fast-assoc run
# ------------------------------------

# Make tassel ready format to run fast association with these PCs
tasselize_me <- subset_nam_FT_PCs[,c(1,3:ncol(subset_nam_FT_PCs))]

# Fix column  names
colnames(tasselize_me)[1] <- "taxa"

# Make column names the first row *tassel format*
tasselize_me[,1] <- as.character(tasselize_me[,1])
tasselize_me <- rbind(colnames(tasselize_me), tasselize_me)

# Create covariate names
tassel_covar_and_data_names <- c("taxa", rep("data", 3), rep("covariate", ncol(tasselize_me)-4))

# Add in covariate, data, taxa IDs
tasselize_me <- rbind(tassel_covar_and_data_names, tasselize_me)

# Export data
write.table(tasselize_me, file = "~/Box Sync/Cornell_PhD/labProjects/hackathons/2019-10-21_hackathon/nam_local_pcs_forTASSEL.txt",
            row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)


# Save image from both of these scripts
save.image(file='nam_cal_local_variance_typeIII.RData')


# -----------------------------
# Plot results from TASSEL run
# -----------------------------

#  Import the data
results <- read.delim("~/Box Sync/Cornell_PhD/labProjects/hackathons/2019-10-21_hackathon/nam_local_global_family_FAresults.txt")

# Remove the 0 pvalues
results <- results[which(results$p > 0),]

# Subset down the file to be compatible with qqman
asi <- results[which(results$Trait == "ASI_BLUP_Sum0607_Buckler2009"),]
dta <- results[which(results$Trait == "Days_To_Anthesis_BLUP_Sum0607_Buckler2009"), ]
dts <- results[which(results$Trait == "Days_To_Silk_BLUP_Sum0607_Buckler2009"), ]

# Plot the results using qqman
library(qqman)

# Make the Manhattan plot on the gwasResults dataset
manhattan(asi, chr="Chr", bp="Pos", snp="Marker", p="p", annotatePval =1.6531e-10,
          main = "y(ASI Buckler 2009) ~ 400k SNPs + 5 gPCs + 660 lPCs")

manhattan(dta, chr="Chr", bp="Pos", snp="Marker", p="p", annotatePval = 1.6531e-10,
          main = "y(DTA Buckler 2009) ~ 400k SNPs + 5 gPCs + 660 lPCs") 

manhattan(dts, chr="Chr", bp="Pos", snp="Marker", p="p", annotatePval = 2.6201e-10,
          main = "y(DTS Buckler 2009) ~ 400k SNPs + 5 gPCs + 660 lPCs")

# Plot QQ plots
qq(asi$p, main = "y(ASI Buckler 2009) ~ 400k SNPs + 5 gPCs + 660 lPCs")
qq(dta$p, main = "y(DTA Buckler 2009) ~ 400k SNPs + 5 gPCs + 660 lPCs")
qq(dts$p, main = "y(DTS Buckler 2009) ~ 400k SNPs + 5 gPCs + 660 lPCs")




