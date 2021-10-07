# ---------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2019-11-04 
#
# Description 
# 	- Compare results of MLM and FA models from TASSEL for
#	  - kinship matrix, global PC, and local PC tests
# ---------------------------------------------------------------

# Load in results
ames2nam_results <- read.delim("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/results/2019_12_04_ames2nam_comparisons/ZeaGBSv27_publicSamples_imputedV5_AGPv4-181023_NAM_only_main200_maf035_FAresults_ames2nam.txt")
nam2nam_results <- read.delim("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/results/2019_12_04_ames2nam_comparisons/ZeaGBSv27_publicSamples_imputedV5_AGPv4-181023_NAM_only_main200_maf035_FAresults_nam2nam.txt")

# Add column for shorter trait names
format_phenotypes <- function(gwas_results){
  gwas_results$short_names <- gsub("_raw.*", "", gwas_results$Trait)
  gwas_results$short_names <- gsub("_Blup.*", "", gwas_results$short_names)
  gwas_results$short_names <- gsub("_Hung.*", "", gwas_results$short_names)
  gwas_results$short_names <- gsub("_BLUP.*", "", gwas_results$short_names)
  gwas_results$short_names <- gsub("_Kump.*", "", gwas_results$short_names)
  gwas_results$short_names_prior <- gsub("southern_leaf_blight", "SLB_resistance", tolower(gwas_results$short_names))
  gwas_results$short_names_prior <- gsub("asi", "flowering_time", gwas_results$short_names_prior)
  gwas_results$short_names_prior <- gsub("days_to_anthesis", "flowering_time", gwas_results$short_names_prior)
  gwas_results$short_names_prior <- gsub("days_to_silk", "flowering_time", gwas_results$short_names_prior)
  gwas_results$short_names_prior <- gsub("chlorophylla", "chlorophyll_a", gwas_results$short_names_prior)
  gwas_results$short_names_prior <- gsub("chlorophyllb", "chlorophyll_b", gwas_results$short_names_prior)
}

# Load prior gene verification list
prior_genes <- read.csv("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/data/verification_loci/verified_genes_traits.csv",
                                   header = T)

# Load packages
library(ggman)
library(ggplot2)
library(GenomicRanges)
library(patchwork)


# ----------------------------
# Check overlapping intervals
# ----------------------------

# Make a function that iterates over all traits and subsets 
holder <- c()
getRanges <- function(df, ld_buffer = 2500000){
  # Store all phenotypes
  traits <- unique(df$trait)
  
  # Remove rows with no start or stop information
  df <- df[!is.na(df$start),]
  
  # Add LD buffer to start and end
  df$start_buffer <- df$start - ld_buffer
  df$stop_buffer <- df$end + ld_buffer

  # Loop over all traits
  for (i in seq(1:length(traits))){
    # Make intervals out of everything
    tmp <- makeGRangesFromDataFrame(df[which(df$trait == traits[i]),], 
                                    start.field = "start_buffer", end.field = "stop_buffer", 
                                    keep.extra.columns = F)

    # Merge overlapping intervals
    tmp <- as.data.frame(reduce(tmp))
    tmp$trait <- rep(traits[i], nrow(tmp))
    
    # Replace negative start sites with 0
    tmp$start[tmp$start<0] <- 0
    
    # add unique intervals to main dataframe with annotation
    holder <- rbind(holder, tmp)
  }
  return(holder)
}

# Use funciton
# Windows do not have consistent size but are close to 5,000,000 bp, with some larger exceptions
prior_ranges <- getRanges(prior_genes)

# Check to see where all GWAS hits are for 1 trait
ggplot(prior_genes[which(prior_genes$trait == "fumarate"),], aes(x = seqid, y = start)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2))


# Make a function that counts the number of overlapping 
#   SNPs between methods
holder <- c()
compare_overlaps <- function(results1, results2, prior_intervals, threshold = 5, maxGap = 100){
  # Get all trait IDs
  all_long <- unique(results1$Trait)
  
  # Iterate through all traits
  for (i in seq(1:length(all_long))){
    # Subset results and phenotypes in both files
    subset_results1 <- results1[which(results1$Trait == all_long[i] & -log10(results1$p) > threshold),]
    subset_results2 <- results2[which(results2$Trait == all_long[i] & -log10(results2$p) > threshold),]

    # Make intervals
    intervals_results1 <- makeGRangesFromDataFrame(subset_results1, start.field = "Pos", end.field = "Pos", keep.extra.columns=F)
    intervals_results2 <- makeGRangesFromDataFrame(subset_results2, start.field = "Pos", end.field = "Pos", keep.extra.columns=F)
    
    # Reduce intervals (combines overlapping ranges)
    intervals_results1 <- reduce(intervals_results1)
    intervals_results2 <- reduce(intervals_results2)
    
    # Get short trait names
    trait_short <- unique(subset_results1$short_names_prior)

    # Subset prior-gene ranges, make GRanges object
    trait_ranges <- prior_intervals[which(prior_intervals$trait == trait_short),]
    trait_ranges <- makeGRangesFromDataFrame(trait_ranges)
    
    # Count the number of overlaps between methods
    overlaps_results1 <- countOverlaps(trait_ranges, intervals_results1, maxgap = maxGap)
    overlaps_results2 <- countOverlaps(trait_ranges, intervals_results2, maxgap = maxGap)

    # Save the total number of windows
    total_windows_results1 <- length(overlaps_results1)
    total_windows_results2 <- length(overlaps_results2)
    
    # Total windows with at least 1 SNP within 100 bp of window
    windows_covered_results1 <- sum(as.logical(overlaps_results1))
    windows_covered_results2 <- sum(as.logical(overlaps_results2))
    
    # Calculate percentage
    percent_results1 <- (windows_covered_results1/total_windows_results1)*100
    percent_results2 <- (windows_covered_results2/total_windows_results2)*100
    
    # Compile results
    together <- cbind(as.character(all_long[i]),
                      windows_covered_results1, total_windows_results1, percent_results1,
                      windows_covered_results2, total_windows_results2, percent_results2)
    
    # Change column names
    colnames(together) <- c("trait", 
                            "windows_captured_ames2nam", "windows_total_ames2nam", "percent_captured_ames2nam",
                            "windows_captured_nam2nam", "windows_total_nam2nam", "percent_captured_nam2nam")
    
    # Add results to one main data.frame
    holder <- rbind(holder, together)
  }
  return(holder)
}

# Use the function
temp <- data.frame(compare_overlaps(ames2nam_results, nam2nam_results, prior_ranges))

# Sort by complexity
temp <- temp[c(6,7,5,9,3,4,8,1,2),]

# Turn into numeric
temp[,2:7] <- data.frame(apply(temp[,2:7], 2, function(x) as.numeric(as.character(x))))


# -------------------
#    Visulization 
# -------------------

#Set directory
setwd("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/results/2019_12_04_ames2nam_comparisons/images")

# Make threshold line
thresh = 6

# Collect all tested phenotypes
traits2plot <- unique(ames2nam_results$short_names)

# For ames2nam and nam2nam phenotypes
for (i in seq(1:length(traits2plot))){
  # ------------------------
  # For ames2nam phenotypes
  # ------------------------
  
  # Shorten name
  results <- ames2nam_results
  
  # Subset trait
  results <- results[which(results$short_names == traits2plot[i]),]
  results <- results[which(results$p >0),]
  
  # Manhattan Plot
  a <- ggman(results, chrom="Chr", bp="Pos", snp="Marker", pvalue="p", pointSize = 1, lineColour = "black",
        title = paste(traits2plot[i],"in NAM \n y ~ x + 3 gPCs + 660 wPCs \n ames2nam"), 
        sigLine = thresh, ymax = 17)
  
  # Setup data
  results$observed <- -log10(results$p)
  results$expected <- -log10(ppoints(length(results$p)))
  
  # QQ Plot
  b <- ggplot(results, aes(x=expected, y=sort(observed,decreasing = T))) +
    geom_point()  +
    xlim(0,5) + ylim(0,17) +
    geom_abline(intercept = 0) +
    ggtitle(paste(traits2plot[i],"in NAM \n y ~ x + 3 gPCs + 660 wPCs \n ames2nam")) +
    labs(x = "Observed", y = "Expected")
  
  # Use patchwork function to plot both
  # print(a+b)
  
  # ----------------------
  # For nam2nam phenotypes
  # ----------------------
  
  # Shorten names
  nam2nam <- nam2nam_results
  
  # Subset trait
  nam2nam <- nam2nam[which(nam2nam$short_names == traits2plot[i]),]
  nam2nam <- nam2nam[which(nam2nam$p >0),]
  
  # Manhattan Plot
  c <- ggman(nam2nam, chrom="Chr", bp="Pos", snp="Marker", pvalue="p", pointSize = 1, lineColour = "black",
             title = paste(traits2plot[i],"in NAM \n y ~ x + 3 gPCs + 660 wPCs \n nam2nam"), 
             sigLine = thresh, ymax = 17)
  
  # Setup data
  nam2nam$observed <- -log10(nam2nam$p)
  nam2nam$expected <- -log10(ppoints(length(nam2nam$p)))
  
  # QQ Plot
  d <- ggplot(nam2nam, aes(x=expected, y=sort(observed,decreasing = T))) +
    geom_point()  +
    xlim(0,5) + ylim(0,17) +
    geom_abline(intercept = 0) +
    ggtitle(paste(traits2plot[i],"in NAM \n y ~ x + 3 gPCs + 660 wPCs \n nam2nam")) +
    labs(x = "Observed", y = "Expected")
  
  # Use patchwork
  all <- (a+b)/(c+d)
  ggsave(paste(i,traits2plot[i], "all_models.png",sep = "_"), plot = all, width = 7.5, height = 9)
}


