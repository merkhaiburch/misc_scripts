# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2021-10-04
# Updated....2022-02-16
#
# Description 
# Calculating pleiotropy across intervals. 
# Run on a lm machine
# ---------------------------------------------------------------

# ----------------------------
# Load in helpful packages
# -----------------------------

library(dplyr)
library(data.table)

# Set global parametes
n_threads <- 60


# ------------------------------------------------------------------------------------------
#                      Gather Data
# ------------------------------------------------------------------------------------------

# -----------------------------
# Gather data from cbsu blfs1
# -----------------------------

# Count Matrices
# mkdir -p /workdir/mbb262/aggregated_counts/aggregated_by_pop_or_trait
# scp mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/mbb262/pleiotropy/interval_data/aggregated_counts/aggregated_by_pop_or_trait/* /workdir/mbb262/aggregated_counts/aggregated_by_pop_or_trait


# Unzip files
# for FILE in *.gz
# do
# bgzip -d --threads 60 ${FILE}
# done


# ---------------------------------------------
# Gather counts of traits within intervals
# counts calculated with Zack's tassel plugin
# ---------------------------------------------

# Different datasets
# - Nam filtered physiological (nfp)
# - Nam permuted physiological (npp)
# - Goodman filtered physiological (gfp)
# - Goodman permuted physiological (gpp)
# - Goodman filtered metabolite (gfm)
# - Goodman permuted metabolite (gpm)
# - Goodman filtered expression (gfe)
# - Goodman permuted expression (gpe)


# ------------------------------------------------------------------------------------------
#                        Non-Permuted data counts 
# ------------------------------------------------------------------------------------------

# Datatbale function to replace NAs with 0's
replace_nas <- function(DT) {
  for (j in seq_len(ncol(DT)))
    set(DT,which(is.na(DT[[j]])),j,0)
  
  return(DT)
}

# NOTE: script 15_aggregate_gwa_counts leaves NAs when there are no GWA hits in any intervals
# we need to replace those NAs with 0's (e.g. there are only 3 sig. snps for a trait, the focal
# interval will have a count of 3 but all other intervals will have NA values) --> ONLY for expression data

# Start loading in data
data_path <- "/workdir/mbb262/aggregated_counts/aggregated_by_pop_or_trait/"

# - Nam filtered physiological (nfp)
nfp <- data.table::fread(paste0(data_path, "nam_filtered_all_traits.txt"))
nfp[is.na(nfp)] <- 0

# - Goodman filtered physiological (gfp)
gfp <- data.table::fread(paste0(data_path, "goodman_filtered_physiological.txt"))
gfp <- gfp %>% select(-contains(".y")) # remove extra Hung phenotypes, they're repeated on accident
colnames(gfp) <- gsub(".x", "", colnames(gfp))
gfp[is.na(gfp)] <- 0

# - Goodman filtered metabolite (gfm)
gfm <- data.table::fread(paste0(data_path, "goodman_filtered_metabolite.txt"))
gfm[is.na(gfm)] <- 0

# - Goodman filtered expression (gfe)
gfe <- data.table::fread(paste0(data_path, "goodman_filtered_expression.txt"), nThread = n_threads)
gfe <- cbind(gfe[,1], replace_nas(gfe[,-1]))
expression_names <- colnames(gfe) # subset out column names from this expression data


# ------------------------------------------------------------------------------------------
#                                       Expression tangent
# ------------------------------------------------------------------------------------------

# For some reason, mapping the same expression data led to different numbers of columns, gather
# column names, intersect them all to get shared names (i.e. same columns) so that pleiotropy 
# is calculated on all the same columns
gather_col_names <- function(file_dir, file_name){
  # Find files with matching string in file name
  find_files <- list.files(file_dir, pattern = file_name, full.names = TRUE)
  
  # Load all into nested list
  perm_list <- lapply(find_files, data.table::fread, nThread = n_threads)
  
  # Gather column names from each list
  gather_names <- lapply(perm_list, colnames)
  
  # remove what permutation part of column names
  for (i in 1:length(gather_names)){
    gather_names[[i]] <- gsub("_permutation_[0-9]", "",gather_names[[i]])
  }
  
  # Intersect all names across lists and return
  df <-Reduce(intersect, gather_names)
  return(df)
}

# Gather goodman expression names shared among all permutations and 
# the biological (non-permuted set) --> length = 116,367
perm_dir <- "/workdir/mbb262/aggregated_counts/aggregated_by_pop_or_trait"
gpe_names <- gather_col_names(file_dir = perm_dir, file_name = "goodman_expression_permutation*")
expression_intersection <- intersect(gpe_names, expression_names)

# number of expression traits we're working with that match across everything:
length(expression_intersection) 


# ------------------------------------------------------------------------------------------
# Calculate pleiotropy metric for non-permuted data
# Metabolite and field traits
# ------------------------------------------------------------------------------------------

# Function to calculate the number of unique traits mapping to each interval ---
calc_pleiotropy_countUnique <- function(aggregated_count_df, data_id_to_paste){

  # Get count across all traits and categories
  # Gather names of traits with associations by SNP
  sub_aggregated_count_df <- aggregated_count_df[,-1]
  traitsMapped <- data.frame(apply(sub_aggregated_count_df, 1, function(x) paste(colnames(sub_aggregated_count_df)[ which(x>0) ], collapse = ",")))
  colnames(traitsMapped) <- "traitsMapped"
  
  # Count of unique traits
  unique_count <- rowSums(sub_aggregated_count_df != 0)
  
  # Rejoin medians with interval IDs
  interval_metrics <- cbind(aggregated_count_df[,1], data.table(unique_count), traitsMapped)

  # Change column names
  colnames(interval_metrics)[2] <- paste0(data_id_to_paste, "_uniqueCount")

  # Return function
  return(interval_metrics)
}

# ------------------------------------------------------------------------------------------
# Calculate pleiotropy metric for non-permuted data: VERSION 2
# Version 2 of the function to count unique traits by trait category
# Field traits
# ------------------------------------------------------------------------------------------

# All populations and trait categories
trait_grouping <- read.csv("/workdir/mbb262/main_phenotypes_melted.csv")

# Gather unique traits to iterate through - NAM
grouping_NAM <- trait_grouping %>% 
  filter(!trait_category %in% c("Metabolite", "Expression") & population == "NAM")
unique_trait_categories_NAM <- unique(grouping_NAM$trait_category)
grouping_NAM$traits <- gsub(";", ".", grouping_NAM$traits) # replace semicolon with period 

# remove expression and metabolite data, filter to only Goodman
grouping_GAP <- trait_grouping %>% 
  filter(!trait_category %in% c("Metabolite", "Expression") & population == "Goodman_Association")
unique_trait_categories_GAP <- unique(grouping_GAP$trait_category)
grouping_GAP$traits <- gsub(";", "_", grouping_GAP$traits) # replace semicolon with period 
grouping_GAP$traits <- paste0(grouping_GAP$traits, "_goodman") # add on "_goodman"

# i <- 1
# aggregated_count_df <- nfp
# pop = "NAM"

calc_pleiotropy_countUnique_fieldTraits <- function(aggregated_count_df, data_id_to_paste, pop){
  # Gather data outside of loop
  outside <- c()
  outside <- data.frame(aggregated_count_df$Interval)
  colnames(outside) <- "Interval"
  
  # select correct population to iterate through
  if (pop == "NAM"){
    unique_trait_categories <- unique_trait_categories_NAM
    grouping <- grouping_NAM
  } else if (pop == "GAP") {
    unique_trait_categories <- unique_trait_categories_GAP
    grouping <- grouping_GAP
  }
  
  # Loop through unique traits and make counts by trait category
  for (i in 1:length(unique_trait_categories)) {
    # Print the trait category
    message(paste0("I am counting trait category: ", unique_trait_categories[i]))
    
    # Select traits within this category in the cross-reference table
    category_sub <- grouping %>% 
      filter(trait_category == unique_trait_categories[i]) %>% 
      select(traits)
    
    # Subset the count df by columns in a particular trait category
    sub_test <- aggregated_count_df %>% select(category_sub$traits)

    # Return column names if trait associates in that interval/SNP
    traitsMapped <- data.frame(apply(sub_test, 1, function(x) paste(colnames(sub_test)[ which(x>0) ], collapse = ",")))

    # Count sum of unique traits in a particular trait category
    unique_count <- rowSums(sub_test != 0)

    # Gather and format data
    together <- cbind(aggregated_count_df[,1], data.frame(unique_count), traitsMapped)
    colnames(together) <- c("Interval", 
                            paste0("unique_count", unique_trait_categories[i]), 
                            paste0("traitsMapped_", unique_trait_categories[i]))

    # Combine outside of loop
    outside <- merge(outside, together, by = "Interval")

  }

  # Get count across all traits and categories
  # Gather names of traits with associations by SNP
  sub_aggregated_count_df <- aggregated_count_df[,-1]
  traitsMapped <- data.frame(apply(sub_aggregated_count_df, 1, function(x) paste(colnames(sub_aggregated_count_df)[ which(x>0) ], collapse = ",")))
  colnames(traitsMapped) <- "traitsMapped"

  # Count of unique traits
  unique_count <- rowSums(sub_aggregated_count_df != 0)

  # Rejoin medians with interval IDs
  interval_metrics <- cbind(aggregated_count_df[,1], data.table(unique_count), traitsMapped)
  
  # Change column names
  colnames(interval_metrics)[2] <- paste0(data_id_to_paste, "_uniqueCount")

  # Merge with individual category counts
  interval_metrics <- merge(interval_metrics, outside, by = "Interval")

  # Return results of function
  return(interval_metrics)
}


# Nam filtered physiological (nfp)
nfp_pleiotropy <- calc_pleiotropy_countUnique_fieldTraits(nfp, "nam_filtered_all", pop = "NAM")

# Goodman filtered physiological (gfp)
gfp_pleiotropy <- calc_pleiotropy_countUnique_fieldTraits(gfp, "goodman_filtered_physiological", pop = "GAP")

# Goodman filtered metabolite (gfm)
gfm_pleiotropy <- calc_pleiotropy_countUnique(gfm, "goodman_filtered_metabolite")

# Goodman filtered expression (gfe) --> CRASHES ON A MM MACHINE, fine on LM machine
dim(gfe)
gfe <- gfe[,..expression_intersection] # subset to shared columns below
dim(gfe) # 116382
gfe_pleiotropy <- calc_pleiotropy_countUnique(gfe, "goodman_filtered_expression")


# Add on identifier for population
colnames(nfp_pleiotropy)[3:19] <- paste0("nfo_", colnames(nfp_pleiotropy)[2:19]) # nam field observed
colnames(gfp_pleiotropy)[3:19] <- paste0("gfo_", colnames(gfp_pleiotropy)[2:19]) # goodman field observed


# ------------------------------------------------------------------------------------------
#       Permuted data counts 
# ------------------------------------------------------------------------------------------

# NOTE: replaced all NAs with 0's in prior aggregate_gwa_counts.R script

# Function to make counts for each permutation
count_nonunique_traits <- function(file_dir, file_name, expression_df = FALSE){
  # Find files with matching string in file name
  find_files <- list.files(file_dir, pattern = file_name, full.names = TRUE)
  
  # Load all into nested list
  perm_list <- lapply(find_files, data.table::fread, nThread = n_threads)
  
  # change column names, subset names to those shared
  if (expression_df == TRUE){
    for (i in 1:length(perm_list)){
      colnames(perm_list[[i]]) <- gsub("_permutation_[0-9]", "", colnames(perm_list[[i]]))
      
      # remove columns that don't have matching names
      perm_list[[i]] <- perm_list[[i]][,..expression_intersection]
    }
  }
  
  # Count number of unique traits within each permutation round
  count_perm <- lapply(perm_list, function(x) cbind(x[,1], rowSums(x[,-1] != 0)))
  print(lapply(perm_list, dim))
  
  # Change colnames names to add permutation info
  for (i in 1:length(find_files)){
    colnames(count_perm[[i]])[2] <- paste0("perm_", i)
  }
  
  # Turn into dataframe
  df <- Reduce(function(x, y) merge(x, y, by = "Interval"), count_perm)

  # Return 
  return(df)
}

# Directory where all files are
perm_dir <- "/workdir/mbb262/aggregated_counts/aggregated_by_pop_or_trait/"

# - Nam permuted physiological (npp)
npp <- count_nonunique_traits(file_dir = perm_dir, 
                    file_name = "nam_physiological_permutation*")

# - Goodman permuted physiological (gpp)
gpp <- count_nonunique_traits(file_dir = perm_dir, 
                    file_name = "goodman_physiological_permutation*")

# - Goodman permuted metabolite (gpm)
gpm <- count_nonunique_traits(file_dir = perm_dir, 
                    file_name = "goodman_metabolite_permutation*")

# - Goodman permuted expression (gpe)
gpe <- count_nonunique_traits(file_dir = perm_dir, 
                    file_name = "goodman_expression_permutation*", expression_df = TRUE)

# NOTE: MIGHT HAVE TO CHANGE

# append on names so merging isn't a nightmare
colnames(npp)[2:11] <- paste0("npp_", colnames(npp)[2:11])
colnames(gpp)[2:11] <- paste0("gpp_", colnames(gpp)[2:11])
colnames(gpm)[2:11] <- paste0("gpm_", colnames(gpm)[2:11])
colnames(gpe)[2:6] <- paste0("gpe_", colnames(gpe)[2:6])


# ------------------------------------------------------------------------------------------
#     Unit tests for permuted and non-permuted data
# ------------------------------------------------------------------------------------------

# Check dimensions between filterd and permuted sets (permuted dimensions printed in above function)
dim(nfp) # nam 75490   122
dim(gfp) # goodman physiological 75490   223
dim(gfm) # goodman metabolite 75490  3874
dim(gfe) # goodman expression --> final # of traits 75490 116382
 
dim(npp) # 75490    11
dim(gpp) # 75490    11
dim(gpm) # 75490    11
dim(gpe) # 75490    6


# They all match in dimensions!!!!! Woo!


# ------------------------------------------------------------------------------------------
#       Combine and export all counts 
# ------------------------------------------------------------------------------------------

# Combine all into single df
all_pleiotropy <- Reduce(function(x, y) merge(x, y, by = "Interval"), 
                         list(nfp_pleiotropy, npp, 
                              gfp_pleiotropy, gpp,
                              gfm_pleiotropy, gpm, 
                              gfe_pleiotropy, gpe))


# Save to file
data.table::fwrite(all_pleiotropy, 
                   file = "/workdir/mbb262/aggregated_counts/aggregated_by_pop_or_trait/goodman_nam_allTraitTypes_gwaCounts_withTraitNames.txt")


# Backing up everything in /workdir/mbb262/aggregated_counts/aggregated_by_pop_or_trait/
# on blfs1 here:
# scp /workdir/mbb262/aggregated_counts/aggregated_by_pop_or_trait/*gz mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/mbb262/results/pleiotropy/interval_data/aggregated_by_pop_or_trait




