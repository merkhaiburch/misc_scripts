# ---------------------------------------
# Merritt Burch
# mbb262@cornell.edu
# 2019-10-10
# Script to do local PCA on the vcf file 
# (feature matrix) on filtered Ames
# SNPs, imputed with Beagle
#
# Data filtering: table sites: MinCount = 2500, MAF 0.0001
# table taxa: MinFreq = 0.55
#
# Script calculates variance components
# ---------------------------------------

# Collector of numbers
all_local_pcs_variance <- c()

# The function
calc_local_variance <- function(genofileID, globalPCs, numLocalPcsOut, windowSize = 1e7){
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
      # Print status
      # cat("\n I am on window ", windows, "/", length(snp_split), "\n")
      
      # Get genotype matrix
      X_list <- snpgdsGetGeno(genofile, snpfirstdim=FALSE, with.id=TRUE, snp.id = unlist(snp_split[windows]), verbose = FALSE)
      
      # Format X matrix
      X <- X_list$genotype
      rownames(X) <- X_list$sample.id
      colnames(X) <- X_list$snp.id
      
      ## to get alternate counts instead reference counts
      X <- 2-X
      
      # Adjust subset of snps using global PCs
      adj_snps <- regress_out(X, globalPCs)
      
      # Calculate and get local PCs
      local_PCs <- prcomp(adj_snps)
      
      # Get the cumulative proportion of variance explained by each PC within a window
      # on a particular chromosome
      tmp <- (local_PCs$sdev^2/sum(local_PCs$sdev^2))[1:numLocalPcsOut]
      collector <- rbind(collector, tmp)
      rownames(collector)[windows] <- paste0("chr_", chrom, "_window_", windows, sep = "")
    }
    # Collect the  PCs in an outside dataframe
    all_local_pcs_variance <- rbind(all_local_pcs_variance, collector)
    colnames(all_local_pcs_variance) <- c("PC1_var", "PC2_var", "PC3_var")
  }
  # Return all the PCs in a big list
  return(all_local_pcs_variance)
}

# Terrible way of getting to the total variance explained by each PC in each window
# Someday need to fix the function to return everything
adj_pcs <- pca_global$PC[,1:5] # Global PCs to adjust by, 5 in this case
local_variance_genome <- calc_local_variance(genofile, adj_pcs, numLocalPcsOut = 3, windowSize = 1e7)

# Get a copy I can work with
local_variance_genome_test <- data.frame(local_variance_genome)

# Melt dataframe into two columns
m_local_var_genome <- reshape::melt(local_variance_genome)

# Rename columns
colnames(m_local_var_genome) <- c("window", "PC", "variance_explained")

# Grep out window number
m_local_var_genome$windowNumber <- gsub(".*window_", "", m_local_var_genome$window)
m_local_var_genome$windowNumber <- as.numeric(m_local_var_genome$windowNumber)


# --------------------------------------
# Plot the variance explained by each PC
# --------------------------------------
 library(ggplot2)

# Break out by facet/chromosome
m_local_var_genome$chr <- gsub("_window_.*", "", m_local_var_genome$window)

# Plot by chromosome and PC
ggplot2::ggplot(m_local_var_genome, aes(x = windowNumber, y = variance_explained, fill = PC)) +
  geom_bar(stat = 'identity')+
  facet_grid(. ~ chr, scales = "free")


# ---------------------------------------
# Investigate SNP desnity across genome
# ---------------------------------------

# why do I only have 1 window across all of chromosome 10?
snp_map_genome <- snpgdsSNPList(genofile)
snp_map_genome_10 <- snp_map_genome[snp_map_genome$chromosome==10,]

# Count number of SNPs per chromosome
tmp <- snp_map_genome
tmp$chromosome <- as.factor(snp_map_genome$chromosome)
temp <- data.frame(table(tmp$chromosome))
colnames(temp) <- c("chrom", "no_snps")

# Plot histograms of snp position across all chromosomes
par(mfrow = c(2,5))
for (chrom in seq(1:10)){
  hist(snp_map_genome[snp_map_genome$chromosome==chrom,]$position,
       main = paste("chromosome", chrom),
       xlab = "SNP Position")
}
par(mfrow = c(1,1))


# ----------------------------------
# Look at the number of SNPs/window
# ----------------------------------

# Get the number of windows and snps per window
# Collector of numbers
snps_per_window <- c()

# The function
snps_per_window_function <- function(genofileID, windowSize = 1e7){
  # Iterate over chromosmes
  for (chrom in seq(1:10)){
    # Get whole SNP list
    snp_map <- snpgdsSNPList(genofileID)
    
    # Subset SNPs by chromosome and window
    windows_genome <- cut_width(snp_map$position[snp_map$chromosome == as.character(chrom)], width = windowSize)
    
    # Subset SNPs by chromosome, window, put them into separate lists to iterate over
    snp_split <- split(x=snp_map$snp.id[snp_map$chromosome == as.character(chrom)], f=windows_genome)
    
    # Another loop to iterate over all windows and calculate local PCs
    collector <- c()
    for (windows in seq(1:length(snp_split))){
      # Get subsetted genotype matrix
      X_list <- snpgdsGetGeno(genofile, snpfirstdim=FALSE, with.id=TRUE, snp.id = unlist(snp_split[windows]), verbose = F)
      
      # Collect SNP facts
      snp_facts <- data.frame(matrix(c(paste0("chr_", chrom, "_window_", windows, sep = ""), chrom, 
                                       length(X_list$snp.id)), nrow = 1, ncol = 3))
      colnames(snp_facts) <- c("window", "chrom", "no_snps")
      
      # Add this data to the collector outside of the loop
      collector <- rbind(collector, snp_facts)
    }
    # Collect the  PCs in an outside dataframe
    snps_per_window <- rbind(snps_per_window, collector)
  }
  
  # Return all the PCs in a big list
  return(snps_per_window)
}

snp_windows <- snps_per_window_function(genofileID = genofile)
snp_windows$no_snps <- as.numeric(as.character(snp_windows$no_snps))

# Plot number of snps per chromosome over all windows
ggplot(snp_windows, aes(x = chrom, y = no_snps)) +
  geom_boxplot() +
  ggtitle("Number of SNPs per window for each chromosome")

# Correlate SNP density (number of snps/window) with variance explained
local_variance_genome_test$sum_variance_3_PCs <- rowSums(local_variance_genome_test[,1:3])
cor(snp_windows$no_snps, local_variance_genome_test$sum_variance_3_PCs)
plot(snp_windows$no_snps, local_variance_genome_test$sum_variance_3_PCs,
     xlab = "Number of SNPs per Window",
     ylab = "Sum of variance explained by 3 local PCs",
     main = "Cumulative variance explained by number of local SNPs/window, R^2 = -0.744")


