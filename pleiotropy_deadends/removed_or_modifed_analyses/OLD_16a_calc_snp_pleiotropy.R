# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-12-19
# Updated....2022-12-19
#
# Description 
# Calculating pleiotropy across intervals at the SNP level
# ---------------------------------------------------------------

# ----------------------------
# Load in helpful packages
# -----------------------------

library(dplyr)
library(data.table)

# Set global parameters
n_threads <- 60


# ------------------------------------------------------------------------------------------
#                      Gather Data
# ------------------------------------------------------------------------------------------

# -----------------------------
# Gather data from cbsu blfs1
# -----------------------------

# Count Matrices
# scp mbb262@cbsublfs1.biohpc.cornell.edu: /workdir/mbb262/aggregated_counts/aggregated_by_pop_or_trait


# unzip Kremling expression results
# bgzip -d --threads 30 


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
data_path <- "/workdir/mbb262/aggregated_counts_snp/aggregated_by_pop_or_trait/"

# - Nam filtered physiological (nfp) --> DONE 2022-12-19 (verified)
nfp <- data.table::fread(paste0(data_path, "nam_filtered_all_traits.txt"))
nfp[is.na(nfp)] <- 0

# - Goodman filtered physiological (gfp) --> DONE 2022-12-19 (verified)
gfp <- data.table::fread(paste0(data_path, "goodman_filtered_physiological.txt"))
gfp <- gfp %>% select(-contains(".y")) # remove extra Hung phenotypes, they're repeated on accident
colnames(gfp) <- gsub(".x", "", colnames(gfp))
gfp[is.na(gfp)] <- 0

# - Goodman filtered metabolite (gfm) --> DONE 2022-12-19 (verified)
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
# ------------------------------------------------------------------------------------------

# Function to calculate the number of unique traits mapping to each intervl
calc_pleiotropy_countUnique <- function(aggregated_count_df, data_id_to_paste){
  
  # Gather names of traits with associations by SNP
  traitsMapped <- data.frame(apply(aggregated_count_df[,-1], 1, function(x) paste(colnames(aggregated_count_df)[ which(x>0) ], collapse = ",")))
  colnames(traitsMapped) <- "traitsMapped"
  
  # Count of unique traits
  unique_count <- rowSums(aggregated_count_df[,-1] != 0)
  
  # Rejoin sum of unique traits with interval IDs
  interval_metrics <- cbind(aggregated_count_df[,1], data.table(unique_count), traitsMapped)
  
  # Change column names
  colnames(interval_metrics)[2] <- paste0(data_id_to_paste, "_uniqueCount")
  
  # Return function
  return(interval_metrics)
}

# Nam filtered physiological (nfp) --> DONE 2022-12-19 (verified)
nfp_pleiotropy <- calc_pleiotropy_countUnique(nfp, "nam_filtered_all")

# Goodman filtered physiological (gfp) --> DONE 2022-12-19 (verified)
gfp_pleiotropy <- calc_pleiotropy_countUnique(gfp, "goodman_filtered_physiological")

# Goodman filtered metabolite (gfm)
gfm_pleiotropy <- calc_pleiotropy_countUnique(gfm, "goodman_filtered_metabolite")

# Goodman filtered expression (gfe) --> CRASHES ON A MM MACHINE, fine on LM machine
dim(gfe)
gfe <- gfe[,..expression_intersection] # subset to shared columns below
dim(gfe) # 116382
gfe_pleiotropy <- calc_pleiotropy_countUnique(gfe, "goodman_filtered_expression")


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


# append on names so merging isn't a nightmare
colnames(npp)[2:11] <- paste0("npp_", colnames(npp)[2:11])
colnames(gpp)[2:11] <- paste0("gpp_", colnames(gpp)[2:11])
colnames(gpm)[2:11] <- paste0("gpm_", colnames(gpm)[2:11])
colnames(gpe)[2:6] <- paste0("gpe_", colnames(gpe)[2:6])


# ------------------------------------------------------------------------------------------
#     Unit tests for permuted and non-permuted data
# ------------------------------------------------------------------------------------------

# Check dimensions between filterd and permuted sets (permuted dimensions printed in above function)
dim(nfp) # nam
dim(gfp) # goodman physiological
dim(gfm) # goodman metabolite
dim(gfe) # goodman expression --> final # of traits 75490 116382

dim(npp)
dim(gpp)
dim(gpm)
dim(gpe)


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
                   file = "/workdir/mbb262/aggregated_counts/aggregated_by_pop_or_trait/goodman_nam_allTraitTypes_gwaCounts.txt")


# Backing up everything in /workdir/mbb262/aggregated_counts/aggregated_by_pop_or_trait/
# on blfs1 here:
# scp /workdir/mbb262/aggregated_counts/aggregated_by_pop_or_trait/*gz mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/mbb262/results/pleiotropy/interval_data/aggregated_by_pop_or_trait




