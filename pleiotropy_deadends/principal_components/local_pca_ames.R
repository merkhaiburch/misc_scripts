# ---------------------------------------
# Merritt Burch
# mbb262@cornell.edu
# 2019-09-25
# Script to do local PCA on the vcf file 
# (feature matrix) on filtered Ames
# SNPs, imputed with Beagle
#
# Data filtering: table sites: MinCount = 2500, MAF 0.0001
# table taxa: MinFreq = 0.55
# ---------------------------------------

# Steps to running local PCA on a genotype (X) matrix
# - Import vcf file
# - get genotype (X) matrix (SNPRelate gives count of reference allele, normally programs giive alterate counts)
# - Scale X matrix with varable number of PCs (1-5)
# - Calculate local PC windows (using list of SNP IDs) on scaled matrix
# - Extract out 1-3 local PCs from each analysis, add to one central dataframe


# 100 cM = 100 Mb

# -------------------------------------------
# Calculate/adjust local PCs with global PCs
# -------------------------------------------

# Load packages
library(SNPRelate)
library(gdsfmt)
library(ggplot2)

# Set directory
setwd("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/ames_tests/Sep25_localPCs")

# Load in Beagle imputed SNPs that were filtered with the 'temp' (aka ames_taxa_gbs_deduplicated.txt, code@bottom)
vcf.fn <- "~/Box Sync/Cornell_PhD/labProjects/hap_gwa/ames_tests/Sep25_localPCs/ZeaGBSv27_20171204_AGPv4_Ames_filteredGenotypes_beagleImputation.vcf.gz"

beagle_vcf <- snpgdsVCF2GDS(vcf.fn, "beagle_filtered.gds", method="copy.num.of.ref", 
                            snpfirstdim=TRUE, ignore.chr.prefix = "S")

# Get summary output
snpgdsSummary("beagle_filtered.gds")

# Check if it was imported correctly
genofile <- snpgdsOpen(beagle_vcf)


# -------------------------------
#     Calculate global PCs
# -------------------------------

# Don't do this, this centers and scales the PCs
# pca <- snpgdsPCA(genofile, num.thread=2) # Centers and scales all in one "hidden" step
# G_scaled <- snpgdsGRM(genofile) ## centered and scaled
# PCA.G_scaled <- decomp(G_scaled$grm, k=5, center = TRUE) # Showing this is the same as 'pca', can cor(pca eigenval, pca.g_scaled eigenvval) = 1

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

# Calculate PCs this way --> centered only
G <- snpgdsGRM(genofile, method="EIGMIX") ## centered only
pca_global <- decomp(G$grm, k=5, center = TRUE)


# --------------------------------
# Function to calculate local PCs
# --------------------------------
# windowsize of 1e7=~300 markers per window

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
all_local_pcs <- c()

# The function
calc_local_pcs <- function(genofileID, globalPCs, numLocalPcsOut, windowSize = 1e7){
  # Iterate over chromosmes
  for (chrom in seq(1:10)){
    
    # Give a status message
    cat("\n I am analyzing chromosome", chrom, "\n")
    
    # Get whole SNP list
    snp_map <- snpgdsSNPList(genofileID)
    
    # Subset SNPs by chromosome and window
    windows_genome <- cut_width(snp_map$position[snp_map$chromosome == as.character(chrom)], width = windowSize)

    # Subset SNPs by chromosome, window, put them into separate lists to iterate over
    snp_split <- split(x=snp_map$snp.id[snp_map$chromosome == as.character(chrom)], f=windows_genome)
    
    # Another loop to iterate over all windows and calculate local PCs
    collector <- c()
    for (windows in seq(1:length(snp_split))){
      #length(snp_split)
      # Print status
      cat("\n I am on window ", windows, "/", length(snp_split), "\n")
      
      # Get genotype matrix
      # Genotype data --> get X matrix
      ## snpfirstdim=FALSE to get a nxp matrix
      ## with.id=TRUE to get sample and SNP information in a list
      X_list <- snpgdsGetGeno(genofileID, snpfirstdim=FALSE, with.id=TRUE, snp.id = unlist(snp_split[windows]))
      
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
      local_PCs <- prcomp(adj_snps)$x
      
      # Save only a subset of the PCs
      subset_local_pcs <- local_PCs[,1:numLocalPcsOut]
      
      # Add local PCs to main dataframe
      collector <- cbind(collector, subset_local_pcs)
      totalCols <- ncol(collector)
      colnames(collector)[totalCols-2] <- paste0("chr_", chrom, "_window_", windows, "_PC1", sep = "")
      colnames(collector)[totalCols-1] <- paste0("chr_", chrom, "_window_", windows, "_PC2", sep = "")
      colnames(collector)[totalCols] <- paste0("chr_", chrom, "_window_", windows, "_PC3", sep = "")
    }
    
    # Collect the  PCs in an outside dataframe
    all_local_pcs <- cbind(all_local_pcs, collector)
  }
  
  # Return all the PCs in a big list
  return(all_local_pcs)
}


# Use function
adj_pcs <- pca_global$PC[,1:5] # Global PCs to adjust by, 5 in this case
localPCs_1Mbp <- calc_local_pcs(genofile, adj_pcs, numLocalPcsOut = 3, windowSize = 1e7)


 # ------------------------------------------------
# Add headers and export in TASSEL-friendly format
# ------------------------------------------------

# Rename for simiplicity
PCs <- as.data.frame(localPCs_1Mbp)

# Add back in the taxa IDs
tasselNames <- cbind(read.gdsn(index.gdsn(genofile, "sample.id")),  PCs)


# -----------------------
#  Add in phenotypic data
# -----------------------

# Load all ames phenotypes and subset down to only what I want
ames_phenos <- read.table("/Users/mbb262-admin/git_projects/haplotype_v_snp_gwas/data/all_Ames_Phenos.txt", header = TRUE)
# ames_phenos <- read.table("C:/Users/merri/git_projects/haplotype_v_snp_gwas/data/all_Ames_Phenos.txt", header = TRUE)
ames_phenos <- ames_phenos[,c(1,6:22)]

# Format PCs for merging
tasselNames$copied_taxa <- gsub(":[0-9]{9}", "", tolower(tasselNames[,1]))
tasselNames$copied_taxa <- gsub("_", "", tasselNames$copied_taxa)

# Merge in 5 global PCs
tasselNames <- cbind(pca_global$PC[,1:5], tasselNames)

# Merge dataframes together
merged_pcs_phenos <- merge(x = ames_phenos, y = tasselNames,
                           by.x = "Entry_ID", by.y = "copied_taxa")

# Rearrage columns
merged_pcs_phenos <- merged_pcs_phenos[ , c(24,2:23,25:ncol(merged_pcs_phenos))]

# Fix column  names
colnames(merged_pcs_phenos)[1] <- "taxa"

# Make column names the first row *tassel format*
merged_pcs_phenos[,1] <- as.character(merged_pcs_phenos[,1])
merged_pcs_phenos <- rbind(colnames(merged_pcs_phenos), merged_pcs_phenos)

# Create covariate names
tassel_covar_and_data_names <- c("taxa", rep("data", ncol(ames_phenos)-1), rep("covariate", ncol(PCs)+5))

# Add in covariate, data, taxa IDs
merged_pcs_phenos <- rbind(tassel_covar_and_data_names, merged_pcs_phenos)

# Export data
write.table(merged_pcs_phenos, file = "ames_local_pcs_forTASSEL.txt",
            row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

# Make smaller version with only romay GDD_DTS data and PCs
romay <- merged_pcs_phenos[,c(1,4,19:ncol(merged_pcs_phenos))]
romay <- na.omit(romay)
write.table(romay, file = "ames_local_pcs_forTASSEL_romay.txt",
            row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)


# --------------------------------
# Look at Type III sums of squares
# https://www.r-bloggers.com/anova-%E2%80%93-type-iiiiii-ss-explained/
# --------------------------------

# Import file with all PCs and GDD_DTS data
dts_pcs <- read.table("local_615_PCs_GDD_DTS_phenos2.txt", header = TRUE) # Need to modify this file (remove tassel headers manually)

# Calculate type III sums of squares
library(car)

# Set options/contrasts that sets sums to zero
options(contrasts = c("contr.sum", "contr.poly"))

# Store the model:
model <- lm(GDD_DTS_BLUP_Romay2013 ~ ., dts_pcs[,-1])

# Finally, call the drop1 function on each model component
# The results give the type III SS, including the p-values from an F-test.
# Not through car package, this is through base R
method1_baseR_typeIII <- drop1(model, .~., test = "F")

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
  ggtitle("Type III Sums of Squares by top 3 local PCs from y(GDD_DTS) ~ 5 gPCs + 615 lPCs")

# Save image from both of these scripts
save.image(file='cal_local_variance_typeIII.RData')


# -----------------------------------------------------
# Calculate Type III sums of squares
# holding out windows with Dong flowering time genes
# -----------------------------------------------------

# Load in Dong 2012 genes from Brandon
dong_genes_with_coordinates <- read.csv("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/ames_tests/Sep25_localPCs/dong_genes_with_coordinates.csv", header = TRUE)

# Get window information
# Function gets coordiantes of window start and stops
get_windows <- function(genofileID, globalPCs, numLocalPcsOut, windowSize = 1e7){
  # Iterate across all chromosomes
  all_chr_formatted <- c()
  for (chrom in seq(1:10)){
    
    # Get whole SNP list
    snp_map <- snpgdsSNPList(genofile)
    
    # Subset SNPs by chromosome and window
    windows_genome <- cut_width(snp_map$position[snp_map$chromosome == as.character(1)], width = 1e7)
    
    # Place in own dataframe
    temp <- data.frame(unique(windows_genome))
    temp <- stringr::str_split_fixed(temp$unique.windows_genome., ",", 2)
    
    # Add in chromosome
    temp <- data.frame(cbind(rep(chrom, nrow(temp)), seq(1:nrow(temp)), temp))
    
    # Remove brackets and other stuff
    temp[,3] <- gsub("\\(", "", temp[,3])
    temp[,3] <- gsub("\\[", "", temp[,3])
    temp[,4] <- gsub("]", "", temp[,4])
    
    # Into numeric
    temp[,3] <- as.numeric(as.character(temp[,3]))
    temp[,4] <- as.numeric(as.character(temp[,4]))

    # Add to collector
    all_chr_formatted <- rbind(all_chr_formatted, temp)
  }
  # Format collector
  all_chr_formatted <- data.frame(all_chr_formatted)
  return(all_chr_formatted)
}

# Use the function
ames_windows <- get_windows(genofile, adj_pcs, numLocalPcsOut = 3, windowSize = 1e7)
colnames(ames_windows) <- c("chr", "window", "window_start", "window_stop")

# Subset SNPs by chromosome, window, put them into separate lists to iterate over
ranges <- merge(ames_windows, dong_genes_with_coordinates,by.x="chr", by.y = "seqid")
intersect_dong_ranges <- ranges[with(ranges, window_start <= start & window_stop >= end),]

# Add column and format names to match my last window notation (to subset)
# Not the best way to do this but I'm out of ideas
intersect_dong_ranges$annotation1 <- paste0("chr_",intersect_dong_ranges$chr,"_window_", intersect_dong_ranges$window, "_PC1")
intersect_dong_ranges$annotation2 <- paste0("chr_",intersect_dong_ranges$chr,"_window_", intersect_dong_ranges$window, "_PC2")
intersect_dong_ranges$annotation3 <- paste0("chr_",intersect_dong_ranges$chr,"_window_", intersect_dong_ranges$window, "_PC3")

# Get only unique windows
remove_from_local1 <- unique(intersect_dong_ranges$annotation1)
remove_from_local2 <- unique(intersect_dong_ranges$annotation2)
remove_from_local3 <- unique(intersect_dong_ranges$annotation3)
remove_from_local <- c(remove_from_local1, remove_from_local2, remove_from_local3)

# Remove those windows from the matrix
local_PCs_DongGenes_removed <- data.frame(localPCs_1Mbp)
myvars <- names(local_PCs_DongGenes_removed) %in% remove_from_local
local_PCs_DongGenes_removed <- local_PCs_DongGenes_removed[!myvars]

# Add back in the taxa IDs
local_PCs_DongGenes_removed <- cbind(read.gdsn(index.gdsn(genofile, "sample.id")),  local_PCs_DongGenes_removed)


# Calculate TypeIII sums of squares

# Load all ames phenotypes and subset down to only what I want
ames_phenos <- read.table("/Users/mbb262-admin/git_projects/haplotype_v_snp_gwas/data/all_Ames_Phenos.txt", header = TRUE)
ames_phenos <- ames_phenos[,c(1,8)]

# Format PCs for merging
local_PCs_DongGenes_removed$copied_taxa <- gsub(":[0-9]{9}", "", tolower(local_PCs_DongGenes_removed[,1]))
local_PCs_DongGenes_removed$copied_taxa <- gsub("_", "", local_PCs_DongGenes_removed$copied_taxa)

# Merge in 5 global PCs
local_PCs_DongGenes_removed <- cbind(pca_global$PC[,1:5], local_PCs_DongGenes_removed)

# Merge dataframes together
merged_pcs_phenos_dong_removed <- merge(x = ames_phenos, y = local_PCs_DongGenes_removed,
                           by.x = "Entry_ID", by.y = "copied_taxa")

# Rearrage columns
merged_pcs_phenos_dong_removed <- merged_pcs_phenos_dong_removed[ , c(1:7,9:ncol(merged_pcs_phenos_dong_removed))]

# Set options/contrasts that sets sums to zero
options(contrasts = c("contr.sum", "contr.poly"))

# Store the model:
model <- lm(GDD_DTS_BLUP_Romay2013 ~ ., merged_pcs_phenos_dong_removed[,-1])

# Finally, call the drop1 function on each model component
# The results give the type III SS, including the p-values from an F-test
method1_baseR_typeIII <- drop1(model, .~., test = "F")

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
  ggtitle("Type III Sums of Squares by top 3 local PCs from y(GDD_DTS) ~ 5 gPCs + 568 lPCs \n 27 windows removed whose locations overlap with Dong 2012 loci")




# --------------------------------------
# Format ames taxa names for merging
#         ONE TIME USE CODE
# --------------------------------------

# Load in taxa ID table from Romay paper (supp table 1)
all_ames_taxa <- read.csv("Romay_2013_supp_table1_Burch_landrace_modified.csv", header = T)

# Format names in second column
all_ames_taxa$copied_taxa <- gsub("-[0-9]$", "", tolower(all_ames_taxa$Accesion.N))
all_ames_taxa <- all_ames_taxa[,c(2,9)] # subset out first two columns

# Format file containing GBS names for easy merging with `all_ames_taxa`
ames_taxa_gbs_names <- read.table("ames_beagle_imputed_taxa_list.txt", header = TRUE)

# Do formatting to filter genotypes
ames_taxa_gbs_names$copy_Taxa <- gsub(":[0-9]{9}", "", tolower(ames_taxa_gbs_names$Taxa)) # Remove gbs names
ames_taxa_gbs_names$copy_Taxa <- gsub("_", "", ames_taxa_gbs_names$copy_Taxa) # Remove underscores
ames_taxa_gbs_names$copy_Taxa <- gsub("\\(", "", ames_taxa_gbs_names$copy_Taxa) # Remove parentheses
ames_taxa_gbs_names$copy_Taxa <- gsub(")", "", ames_taxa_gbs_names$copy_Taxa) # Remove parentheses
ames_taxa_gbs_names$copy_Taxa <- gsub("-[0-9]$", "", ames_taxa_gbs_names$copy_Taxa) # For some PI accessions with -1 thru -9

# Merge two datasets together
# sanity check, it's n = 2813 not 2815 (Cinta verified 2815)
names_ames <- merge(x = all_ames_taxa, y = ames_taxa_gbs_names, by.x = "copied_taxa", by.y = "copy_Taxa", all = T)
length(unique(names_ames$copied_taxa))
length(unique(names_ames$Accesion.N))

# Subset names file based on unique names in copied_taxa column and remove missing values
# Did Accesion.N column before, removed too many taxa for my purpose
Taxa <- names_ames[!duplicated(names_ames$copied_taxa),]
colnames(Taxa) <- c("easy_merge_name", "Romay_Accesion.N", "Taxa")
temp <- data.frame(Taxa$Taxa)
colnames(temp) <- "Taxa"
temp <- na.omit(temp) # n=2944

# Export new taxa list
write.table(temp, "ames_taxa_gbs_deduplicated.txt", row.names = F, quote = F)


