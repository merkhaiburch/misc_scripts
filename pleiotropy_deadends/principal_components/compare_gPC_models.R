# ---------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2020-01-06 
#
# Description 
# 	- Compare results of FA models from TASSEL for
#	  - number of global PCs, and method of calculation
#   - Test the impact of founder SNPs and all NAM SNPs using shared
#   - set of 40k SNPs calculated through transfer from ames (ames2nam)
#	  - in fixed and gene based windows

# Ames to NAM PC transferred models 
#   5) y ~ SNP + 3 gPCs + (all NAM, fixed window size)
#   6) y ~ SNP + 3 gPCs + (all NAM, 181 genes/window)
#   7) y ~ SNP + 26 gPCs + (all NAM, fixed window size)
#   8) y ~ SNP + 26 gPCs + (all NAM, 181 genes/window)

# NAM to NAM PC transferred models
# Do the founders alone create good covariates/PCs?
#   9) y ~ SNP + 3 gPCs + (NAM founders only, fixed window size)
#   10) y ~ SNP + 3 gPCs + (NAM founders only, 181 genes/window)
#   11) y ~ SNP + 26 gPCs + (NAM founders only, fixed window size)
#   12) y ~ SNP + 26 gPCs + (NAM founders only, 181 genes/window)

# Do we need all NAM families create good covariates/PCs?
#   13) y ~ SNP + 3 gPCs + (all NAM, fixed window size)
#   14) y ~ SNP + 3 gPCs + (all NAM, 181 genes/window)
#   15) y ~ SNP + 26 gPCs + (all NAM, fixed window size)
#   16) y ~ SNP + 26 gPCs + (all NAM, 181 genes/window)
# ---------------------------------------------------------------


# Format model names (used further down)
model05 = "y~SNP+3gPCs+(ames to nam, all NAM, fixed window)"
model06 = "y~SNP+3gPCs+(ames to nam, all NAM, 181genes/window)"
model07 = "y~SNP+26gPCs+(ames to nam, all NAM, fixed window)"
model08 = "y~SNP+26gPCs+(ames to nam, all NAM, 181genes/window)"

# NAM to NAM PC transferred models
# Do the founders alone creategood covariates/PCs?
model09 = "y~SNP+3gPCs+(NAM founders to all NAM, fixed window)"
model10 =  "y~SNP+3gPCs+(NAM founders to all NAM, 181genes/window)"
model11 =  "y~SNP+26gPCs+(NAM founders to all NAM, fixed window)"
model12 =  "y~SNP+26gPCs+(NAM founders to all NAM, 181genes/window)"

# Do we need all NAM families create good covariates/PCs?
model13 =  "y~SNP+3gPCs+(all NAM to all NAM, fixed window)"
model14 =  "y~SNP+3gPCs+(all NAM to all NAM, 181genes/window)"
model15 =  "y~SNP+26gPCs+(all NAM to all NAM, fixed window)"
model16 = "y~SNP+26gPCs+(all NAM to all NAM, 181genes/window)"

model_names <- c(model05, model06, model07, model08, model09, model10, model11, model12, model13, model14, model15, model16)

# Set directory
setwd("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/results/2020_01_03_ames2nam_comparisons")

# Load packages
library(ggman)
library(ggplot2)
library(GenomicRanges)
library(patchwork)
library(dplyr)

# Load in all files to the same list
filenames <- list.files(pattern="*.txt", full.names=TRUE)
results <- lapply(filenames, read.delim)
names(results) <- substr(filenames, 3, nchar(filenames)-14)

# Bug testing, looking at one dataframe alone
temp <- read.delim("model6_ames2nam_3gPCs_allNAM_geneWindow_FAresults_debug0110.txt")

setwd("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/results/2020_01_30_ames2nam_newMarkers")
model <- read.delim("manual.txt")

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
  return(gwas_results)
}

# Use lapply
fa_results <- lapply(results, format_phenotypes)
temp_formatted <- format_phenotypes(model)

# Load prior gene verification list
prior_genes <- read.csv("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/data/verification_loci/verified_genes_traits.csv",
                        header = T)


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

# Make a function that counts the number of overlapping 
#   SNPs between methods
holder <- c()
compare_overlaps <- function(results1, prior_intervals = prior_ranges, threshold = 5, maxGap = 100){
  # Get all trait IDs
  all_long <- unique(results1$Trait)

  # Iterate through all traits
  for (i in seq(1:length(all_long))){
    # Subset results and phenotypes in both files
    subset_results1 <- results1[which(results1$Trait == all_long[i] & -log10(results1$p) > threshold),]
    print("here")
    # Make intervals
    intervals_results1 <- makeGRangesFromDataFrame(subset_results1, start.field = "Pos", end.field = "Pos", keep.extra.columns=F)
    print("here 2")
    # Reduce intervals (combines overlapping ranges)
    intervals_results1 <- reduce(intervals_results1)
    
    # Get short trait names
    trait_short <- unique(subset_results1$short_names_prior)
    
    # Subset prior-gene ranges, make GRanges object
    trait_ranges <- prior_intervals[which(prior_intervals$trait == trait_short),]
    trait_ranges <- makeGRangesFromDataFrame(trait_ranges)
    
    # Count the number of overlaps between methods
    overlaps_results1 <- countOverlaps(trait_ranges, intervals_results1, maxgap = maxGap)
    
    # Save the total number of windows
    total_windows_results1 <- length(overlaps_results1)
    
    # Total windows with at least 1 SNP within 100 bp of window
    windows_covered_results1 <- sum(as.logical(overlaps_results1))
    
    # Calculate percentage
    percent_results1 <- (windows_covered_results1/total_windows_results1)*100
    
    # Compile results (full version)
    # together <- cbind(as.character(all_long[i]),
    #                   windows_covered_results1, total_windows_results1, percent_results1)
    # Compile results (reduced version)
    together <- cbind(as.character(all_long[i]), percent_results1)
    
    # Change column names (full version)
    # colnames(together) <- c("trait", paste("windows_captured", id, sep = "_"), 
    #                         paste("windows_total", id, sep = "_"), 
    #                         paste("percent_captured", id, sep = "_"))
    # Reduced version
    colnames(together) <- c("trait", "percent_captured")
    
    # Add results to one main data.frame
    holder <- rbind(holder, together)
  }
  return(t(holder))
}

# Use the function
all <- lapply(fa_results, compare_overlaps)
all <- compare_overlaps(temp_formatted)

# Combine all lists into the same dataframe
result_percentages <- data.frame(do.call(rbind, all))

# Make first row colnames, remove every other row, not elegant but it works
colnames(result_percentages) <- unlist(result_percentages[1,])
result_percentages <- result_percentages[-c(1,3,5,7,9,11,13,15,17,19,21,23),]
result_percentages <- cbind(as.factor(model_names), result_percentages)

# Turn into numeric
result_percentages[,2:9] <- data.frame(apply(result_percentages[,2:9], 2, function(x) as.numeric(as.character(x))))

write.table(result_percentages, 
            file = "~/Box Sync/Cornell_PhD/labProjects/hap_gwa/results/2020_01_03_ames2nam_comparisons/percentages.txt",
            quote = T, row.names = F)

# -------------------
#    Visulization 
# -------------------

#Set directory
setwd("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/results/2020_01_03_ames2nam_comparisons/images")

# Make threshold line
thresh = 6

# Collect all tested phenotypes
traits2plot <- unique(fa_results[[1]]$short_names)
traits2plot <- unique(temp_formatted$short_names)
# For ames2nam and nam2nam phenotypes
plot_results <- function(df){
  for (i in seq(1:length(traits2plot))){

    # Subset by trait, remove all 0 p values
    gwa_results <- df[which(df$short_names == traits2plot[i]),]
    gwa_results <- gwa_results[which(gwa_results$p >0),]
    
    # # Get ID of model, multiple dataframes in a list
    # for(m in 1:seq_along(model_names)){
    #   id <- model_names[m]
    # }
    # plot_name <- paste(traits2plot[i]," in NAM \n ", id, sep = "")
    
    # One data.frame
    plot_name <- paste(traits2plot[i]," in NAM \n " ,sep = "")
    
    # Manhattan Plot
    source('~/git_projects/haplotype_v_snp_gwas/src/R/random_scripts/manhattan_plot.R')
    a <- manhattan_plot(gwa_results, sig_threshold = 7, ylim = 17, title = plot_name)

    # Setup data
    gwa_results$observed <- -log10(gwa_results$p)
    gwa_results$expected <- -log10(ppoints(length(gwa_results$p)))
    
    # QQ Plot
    b <- ggplot(gwa_results, aes(x=expected, y=sort(observed,decreasing = T))) +
      geom_point(size = 1)  +
      xlim(0 , 5) + 
      geom_abline(intercept = 0) +
      scale_y_continuous(expand = c(0,0), limits = c(0, 17)) +
      scale_size_continuous(range = c(0,3)) +
      labs(x = "Observed", y = "Expected") +
      theme_minimal() +
      theme( 
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(size = 9, vjust = 0.5),
        axis.text.y = element_text(size = 9, vjust = 0.5))
    
    # Use patchwork
    all <- (a+b)
    ggsave(paste(i,traits2plot[i], "FA_models.png",sep = "_"), plot = all, width = 7.5, height = 4)
  }
}

# Use function to iterate through list with many dfs
for(m in seq_along(fa_results)){
  # Say what model we're on
  print(names(fa_results[m]))
  model_id <- names(fa_results[m])
  # Loop over all traits
  for (i in seq(1:length(traits2plot))){
    
    # Subset by trait, remove all 0 p values
    gwa_results <- do.call(rbind.data.frame, fa_results[m])
    gwa_results <- gwa_results[which(gwa_results$short_names == traits2plot[i]),]
    gwa_results <- gwa_results[which(gwa_results$p >0),]

    # Get ID of model
    plot_name <- paste(traits2plot[i]," in NAM \n ", model_id, sep = "")
    print(plot_name)

    # Manhattan Plot
    source('~/git_projects/haplotype_v_snp_gwas/src/R/random_scripts/manhattan_plot.R')
    a <- manhattan_plot(gwa_results, sig_threshold = 7, ylim = 17, title = plot_name)

    # Setup data
    gwa_results$observed <- -log10(gwa_results$p)
    gwa_results$expected <- -log10(ppoints(length(gwa_results$p)))

    # QQ Plot
    b <- ggplot(gwa_results, aes(x=expected, y=sort(observed,decreasing = T))) +
      geom_point(size = 1)  +
      xlim(0 , 5) +
      geom_abline(intercept = 0) +
      scale_y_continuous(expand = c(0,0), limits = c(0, 17)) +
      scale_size_continuous(range = c(0,3)) +
      labs(x = "Observed", y = "Expected") +
      theme_minimal() +
      theme(
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(size = 9, vjust = 0.5),
        axis.text.y = element_text(size = 9, vjust = 0.5))

    # Use patchwork
    all <- (a+b)
    ggsave(paste(m,i,plot_name, "FA_models.png", sep = "_"), plot = all, width = 7.5, height = 4)
  }
}

plot_results(temp_formatted)                           
                                
# Subset by trait, remove all 0 p values
df = model
gwa_results <- df[which(df$Trait == "Days_To_Anthesis_BLUP_Sum0607_Buckler2009"),]
gwa_results <- df[which(df$Trait == "Days_To_Silk_BLUP_Sum0607_Buckler2009"),]
# gwa_results <- gwa_results[which(gwa_results$p >0),]

# One data.frame
# plot_name <- paste(traits2plot[i]," in NAM \n " ,sep = "")

# Manhattan Plot
# source('~/git_projects/haplotype_v_snp_gwas/src/R/random_scripts/manhattan_plot.R')
a <- manhattan_plot(gwa_results, sig_threshold = 7, ylim = 17, title = "DTS")

# Setup data
gwa_results$observed <- -log10(gwa_results$p)
gwa_results$expected <- -log10(ppoints(length(gwa_results$p)))

# QQ Plot
b <- ggplot(gwa_results, aes(x=expected, y=sort(observed,decreasing = T))) +
  geom_point(size = 1)  +
  xlim(0 , 5) + 
  geom_abline(intercept = 0) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 17)) +
  scale_size_continuous(range = c(0,3)) +
  labs(x = "Observed", y = "Expected") +
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 9, vjust = 0.5),
    axis.text.y = element_text(size = 9, vjust = 0.5))

# Use patchwork
(a+b)

                                
                                