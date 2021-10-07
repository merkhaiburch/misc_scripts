# ---------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2020-02-21 
#
# Description 
#   - Script to count the number of overlaps between your GWAS
#     results and prior verified gene intervals and return a table.
# ---------------------------------------------------------------

# Load Packages
library(GenomicRanges)


# Formats long phenotype names into short names
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


# Creates a Genomic Ranges object that contains prior known genes
# prior_verified_gene_intervals <- function(gene_ld_buffer, gene_type){
#   
#   if (gene_type == "prior"){
#     # Load prior gene verification list
#     # genes <- read.csv("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/data/verification_loci/verified_genes_traits.csv",
#     #                         # header = T)
#     genes <- read.csv("/workdir/mbb262/verified_genes_traits.csv", header = T)
#     
#     # Function takes prior known genes involved in differnt traits
#     # creates a genomic ranges data frame, adds a ld_buffer to the 
#     # start and end of the gene, and merges overlapping ranges
#     holder <- c()
#     getRanges <- function(df = genes, gene_ld_buffer){
#       
#       # Store all phenotypes
#       traits <- unique(df$trait)
#       
#       # Remove rows with no start or stop information
#       df <- df[!is.na(df$start),]
#       
#       # Add LD buffer to start and end of each gene
#       df$start_buffer <- df$start - gene_ld_buffer
#       df$stop_buffer <- df$end + gene_ld_buffer
#       
#       # Loop over all traits in prior gene list
#       for (i in seq(1:length(traits))){
#         # Make intervals out of everything
#         tmp <- makeGRangesFromDataFrame(df[which(df$trait == traits[i]),], 
#                                         start.field = "start_buffer", end.field = "stop_buffer", 
#                                         keep.extra.columns = F)
#         
#         # Merge overlapping intervals
#         tmp <- as.data.frame(reduce(tmp))
#         tmp$trait <- rep(traits[i], nrow(tmp))
#         
#         # Replace negative start sites with 0
#         tmp$start[tmp$start<0] <- 0
#         
#         # add unique intervals to main dataframe with annotation
#         holder <- rbind(holder, tmp)
#         }
#       return(holder)
#     }
#   }
#     
#   if (gene_type == "random"){
#     # Random set of genes
#     genes <- read.csv("/workdir/mbb262/random_genes_equal_length_priors.csv", header = T)
#     # Function takes prior known genes involved in differnt traits
#     # creates a genomic ranges data frame, adds a ld_buffer to the 
#     # start and end of the gene, and merges overlapping ranges
#     holder <- c()
#     getRanges <- function(df = genes, gene_ld_buffer){
# 
#       # Add LD buffer to start and end of each gene
#       df$start_buffer <- df$start - gene_ld_buffer
#       df$stop_buffer <- df$end + gene_ld_buffer
# 
#       # Make GRanges object from all genes
#       tmp <- makeGRangesFromDataFrame(df, start.field = "start_buffer", end.field = "stop_buffer", 
#                                         keep.extra.columns = F)
# 
#       # Merge overlapping intervals
#       tmp <- as.data.frame(reduce(tmp))
#       
#       # Replace negative start sites with 0
#       tmp$start[tmp$start<0] <- 0
#         
#       # add unique intervals to main dataframe with annotation
#       holder <- rbind(holder, tmp)
#       return(holder)
#       }
#     }
# 
#   # Use and return function
#   gene_intervals <- getRanges(genes, gene_ld_buffer)
#   return(gene_intervals)
# }


# Counts number of significant overlaps between GWAS results and the prior known gene list
# but using p-value quantiles of 0.1 and 1% tails of the distribution (instead of absolute p-values)
count_sig_gwas_overlaps_quantile <- function(gwas_results, prior_ranges, max_snp_overlap, df_id){
  # Make a function that counts the number of overlapping 
  #   SNPs between methods
  holder_overlaps <- c()
  compare_overlaps <- function(results = gwas_results, prior_intervals = prior_ranges, overlap = max_snp_overlap, id = df_id){
    # Get all long trait names from Trait columm from Fast Association
    all_long <- unique(results$Trait)
    
    # Iterate through all traits in results
    for (i in seq(1:length(all_long))){
      
      # Subset GWAS results by phenotype (loop iterator)
      subset_results <- results[which(results$Trait == all_long[i]),]
      
      # Calculate two significance thresholds to subset by
      sig_threshold_1 <- quantile(subset_results$p, probs = 0.01)
      sig_threshold_0.1 <- quantile(subset_results$p, probs = 0.0001)
      
      # Print cutoffs from before
      print(all_long[i])
      print("Percentiles")
      print(-log10(quantile(subset_results$p, probs = c(0.01, 0.0001))))

      # Subset data by pvalue
      subset_results_1 <- subset_results[which(subset_results$p <= sig_threshold_1),]
      subset_results_0.1 <- subset_results[which(subset_results$p <= sig_threshold_0.1),]
      
      # Format subsetted results
      subset_results_1 <- format_phenotypes(subset_results_1)
      
      # Subset prior-gene ranges to only the phenotype in question (loop iterator), 
      # make GRanges object
      # Get short trait names (from format_phenotypes step)
      trait_short <- unique(subset_results_1$short_names_prior)
      trait_ranges <- prior_intervals[which(prior_intervals$trait == trait_short),]
      trait_ranges <- makeGRangesFromDataFrame(trait_ranges)

      # If subset yields no significant results for 1% (0.01), quantile paste "NS" into row
      if (nrow(subset_results_1) == 0){
        total_windows_results_1 <- length(trait_ranges)
        windows_covered_results_1 <- "NA"
        percent_results_1 <- 0
      } else {

        # Make intervals (Genomic Ranges) object that is the length of SNP
        intervals_results_1 <- makeGRangesFromDataFrame(subset_results_1, start.field = "Pos", 
                                                      end.field = "Pos", keep.extra.columns=F)
        
        # Count the number of overlaps between methods
        overlaps_results_1 <- countOverlaps(trait_ranges, intervals_results_1, maxgap = overlap)
        
        # Total windows with at least 1 SNP within 100 bp of window
        windows_covered_results_1 <- sum(as.logical(overlaps_results_1))
        
        # Save the total number of windows
        total_windows_results_1 <- length(overlaps_results_1)
        
        # Calculate percentage of overlapping windows
        percent_results_1 <- (windows_covered_results_1/total_windows_results_1)*100
      }

      # If subset yields no significant results for 0.1% (0.001), quantile paste "NS" into row
      if (nrow(subset_results_0.1) == 0){
        total_windows_results_0.1 <- length(trait_ranges)
        windows_covered_results_0.1 <- "NA"
        percent_results_0.1 <- 0
      } else {
        
        # Make intervals (Genomic Ranges) object that is the length of SNP
        intervals_results_0.1 <- makeGRangesFromDataFrame(subset_results_0.1, start.field = "Pos", 
                                                      end.field = "Pos", keep.extra.columns=F)
        
        # Count the number of overlaps between methods
        overlaps_results_0.1 <- countOverlaps(trait_ranges, intervals_results_0.1, maxgap = overlap)
        
        # Total windows with at least 1 SNP within 100 bp of window
        windows_covered_results_0.1 <- sum(as.logical(overlaps_results_0.1))
        
        # Save the total number of windows
        total_windows_results_0.1 <- length(overlaps_results_0.1)
        
        # Calculate percentage of overlapping windows
        percent_results_0.1 <- (windows_covered_results_0.1/total_windows_results_0.1)*100
      }
      
      # Compile results (full version)
      together <- cbind(as.character(all_long[i]),
                        as.numeric(as.character(windows_covered_results_1)), 
                        as.numeric(as.character(percent_results_1)),
                        as.numeric(as.character(windows_covered_results_0.1)), 
                        as.numeric(as.character(percent_results_0.1)),
                        as.numeric(as.character(total_windows_results_1)))
      colnames(together) <- c("trait", 
                              paste0("genes_captured_1_", id),
                              paste0("percent_genes_captured_1_", id),
                              paste0("genes_captured_0.1_", id),
                              paste0("percent_genes_captured_0.1_", id),
                              paste0("genes_total_", id))

      
      # Add results to one main data.frame
      holder_overlaps <- rbind(holder_overlaps, together)
    }
    # Return Function results
    return(data.frame(holder_overlaps))
    # End of function
  }
  # Use the above function
  count <- compare_overlaps(results = gwas_results, prior_intervals = prior_ranges,
                            overlap = max_snp_overlap, id = df_id)
  return(count)
}


# Creates a Genomic Ranges object that contains prior known genes
prior_verified_gene_intervals <- function(gene_ld_buffer, gene_type){
  
  if (gene_type == "random"){
    # Random set of genes
    genes <- read.csv("/workdir/mbb262/random_genes_equal_length_priors.csv", header = T)
  }
  
  if (gene_type == "prior"){
    # genes <- read.csv("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/data/verification_loci/verified_genes_traits.csv",
    #                         # header = T)
    genes <- read.csv("/workdir/mbb262/verified_genes_traits.csv", header = T)
  }
  
  # Function takes prior known genes involved in differnt traits
  # creates a genomic ranges data frame, adds a ld_buffer to the 
  # start and end of the gene, and merges overlapping ranges
  holder <- c()
  getRanges <- function(df = genes, gene_ld_buffer){
    
    # Store all phenotypes
    traits <- unique(df$trait)
    
    # Remove rows with no start or stop information
    df <- df[!is.na(df$start),]
    
    # Add LD buffer to start and end of each gene
    df$start_buffer <- df$start - gene_ld_buffer
    df$stop_buffer <- df$end + gene_ld_buffer
    
    # Loop over all traits in prior gene list
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
  # Use and return function
  gene_intervals <- getRanges(genes, gene_ld_buffer)
  return(gene_intervals)
}


# ----------------------
# Old unused function(s)
# ----------------------
# Counts number of significant overlaps between GWAS results and the prior known gene list
count_sig_gwas_overlaps <- function(gwas_results, prior_ranges, log10_sig_threshold, max_snp_overlap, df_id){
  # Make a function that counts the number of overlapping 
  #   SNPs between methods
  holder_overlaps <- c()
  compare_overlaps <- function(results = gwas_results, prior_intervals = prior_ranges, sig_threshold = log10_sig_threshold, overlap = max_snp_overlap, id = df_id){
    # Get all long trait names from Trait columm from Fast Association
    all_long <- unique(results$Trait)
    
    # Iterate through all traits in results
    for (i in seq(1:length(all_long))){
      
      # Subset GWAS results by phenotype (loop iterator) and 0 p-values
      subset_results <- results[which(results$Trait == all_long[i]),]
      
      # Get short trait names (from format_phenotypes step)
      trait_short <- unique(subset_results$short_names_prior)
      
      # Subset prior-gene ranges to only the phenotype in question (loop iterator), 
      # make GRanges object
      trait_ranges <- prior_intervals[which(prior_intervals$trait == trait_short),]
      trait_ranges <- makeGRangesFromDataFrame(trait_ranges)
      
      # Subset data by pvalue
      subset_results <- subset_results[which(-log10(subset_results$p) > sig_threshold),]
      
      # If subset yields no significant results, paste "NS" into row
      if (nrow(subset_results) == 0){
        total_windows_results <- length(trait_ranges)
        windows_covered_results <- "NA"
        percent_results <- 0
      } else {
        
        # Make intervals (Genomic Ranges) object that is the length of SNP
        intervals_results <- makeGRangesFromDataFrame(subset_results, start.field = "Pos", 
                                                      end.field = "Pos", keep.extra.columns=F)
        
        # Count the number of overlaps between methods
        overlaps_results <- countOverlaps(trait_ranges, intervals_results, maxgap = overlap)
        
        # Total windows with at least 1 SNP within 100 bp of window
        windows_covered_results <- sum(as.logical(overlaps_results))
        
        # Save the total number of windows
        total_windows_results <- length(overlaps_results)
        
        # Calculate percentage of overlapping windows
        percent_results <- (windows_covered_results/total_windows_results)*100
      }
      
      # Compile results (full version)
      together <- cbind(as.character(all_long[i]),
                        as.numeric(as.character(windows_covered_results)), 
                        as.numeric(as.character(total_windows_results)), 
                        as.numeric(as.character(percent_results)))
      colnames(together) <- c("trait", paste("windows_captured", id, sep = "_"),
                              paste("windows_total", id, sep = "_"),
                              paste("percent_captured", id, sep = "_"))
      
      # Compile results (reduced version)
      # together <- cbind(as.character(all_long[i]), percent_results)
      # colnames(together) <- c("trait", "percent_captured")
      
      # Add results to one main data.frame
      holder_overlaps <- rbind(holder_overlaps, together)
    }
    # Return Function results
    return(data.frame(holder_overlaps))
    # End of function
  }
  # Use the above function
  count <- compare_overlaps(results = gwas_results, prior_intervals = priors,
                            sig_threshold = log10_sig_threshold, overlap = max_snp_overlap, id = df_id)
  return(count)
}
