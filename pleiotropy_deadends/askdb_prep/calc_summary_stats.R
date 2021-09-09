# ---------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2020-07-08
#
# Description 
#   - Calculate summary statistics for NAM & 282 data for askDB
# ---------------------------------------------------------------

# Load in packages
library(dplyr)
library(tidyverse)

# Load in files
metabolites_df <- read.csv("~/git_projects/haplotype_v_snp_gwas/data/zhou_metabolites_Assoc282_Phenos.csv", header = TRUE)
goodman_df <- read.csv("~/git_projects/haplotype_v_snp_gwas/data/all_Assoc282_Phenos.csv", header = TRUE)
nam_df <- read.csv("~/git_projects/haplotype_v_snp_gwas/data/all_NAM_phenos.csv", header = TRUE)

# Remove annotation columns
metabolites_df <- metabolites_df[,-c(1:11)]
goodman_df <- goodman_df[,-c(1:11)]
nam_df <- nam_df[,-c(1:8)]


# Calculate summary statistics for metabolites
stats <- lapply(1:ncol(metabolites_df), function(i) {
  x <- metabolites_df[[i]]
  x <- x[!is.na(x)]

  summary_df <- tibble::tibble(
    trait  = names(metabolites_df)[i],
    n      = length(x),
    mean   = mean(x),
    median = median(x),
    sd     = sd(x),
    min    = min(x),
    max    = max(x)
  )

  return(summary_df)
})
stats <- do.call("rbind", stats)


# Calculate summary statistics for goodman data
stats_goodman <- lapply(1:ncol(goodman_df), function(i) {
  x <- goodman_df[[i]]
  x <- x[!is.na(x)]

  summary_df <- tibble::tibble(
    trait  = names(goodman_df)[i],
    n      = length(x),
    mean   = mean(x),
    median = median(x),
    sd     = sd(x),
    min    = min(x),
    max    = max(x)
  )

  return(summary_df)
})
stats_goodman <- do.call("rbind", stats_goodman)

# Calculate summary statistics for NAM data
stats_nam <- lapply(1:ncol(nam_df), function(i) {
  x <- nam_df[[i]]
  x <- x[!is.na(x)]
  
  summary_df <- tibble::tibble(
    trait  = names(nam_df)[i],
    n      = length(x),
    mean   = mean(x),
    median = median(x),
    sd     = sd(x),
    min    = min(x),
    max    = max(x)
  )
  
  return(summary_df)
})
stats_nam <- do.call("rbind", stats_nam)

# Export both files as a csv
write.csv(stats, "goodman_metabolite_pheno_stats.csv", row.names = FALSE, quote = FALSE)
write.csv(stats_goodman, "~/git_projects/haplotype_v_snp_gwas/data/metadata/phenostat_physiological_goodman.csv", 
          row.names = FALSE, quote = FALSE)
write.csv(stats_nam, "~/git_projects/haplotype_v_snp_gwas/data/metadata/phenostat_nam.csv", 
          row.names = FALSE, quote = FALSE)


# ---------------------------------------------
# Create summary statistics for Kremling data
# ---------------------------------------------

# Load in files
raw_files <- list.files("/workdir/mbb262/kremling_raw_count_v4_hapmap321taxaid/", pattern = ".csv", full.names = TRUE)
raw_files_2 <- list.files("/workdir/mbb262/kremling_raw_count_v4_hapmap321taxaid/", pattern = ".csv")
pure_names <- gsub("_kremling_formatted_v4_hapmapids", "", raw_files_2)
pure_names <- gsub(".csv", "", pure_names)

# Calculate summary statistics for all other data
stats <- lapply(1:ncol(raw_phenos), function(i) {
  x <- raw_phenos[[i]]
  x <- x[!is.na(x)]

  summary_df <- tibble::tibble(
    trait  = names(raw_phenos)[i],
    n      = length(x),
    mean   = mean(x),
    median = median(x),
    sd     = sd(x),
    min    = min(x),
    max    = max(x)
  )

  return(summary_df)
})

# Loop through all 7 tissue files
for (i in seq_len(length(raw_files))){
  # Read in file
  raw_phenos <- data.table::fread(
        input = raw_files[i]
    )

  # Remove annotation columns
  raw_phenos <- raw_phenos[,-1]

  # Run summary stat function
  stat_results <- do.call("rbind", stats)

  # Paste tissue name in front of gene ID
  stat_results$trait <- paste0(pure_names[i], "_", stat_results$trait)

  # Export both files as a csv
  name <- paste0("pheno_stats_", raw_files_2[i])
  write.csv(stat_results, file = name, row.names = FALSE, quote = FALSE)
}


# ---------------------------------
# Make metadata csv for all traits
# ---------------------------------

metadata <- read.csv("/workdir/mbb262/kremling_metadata.csv", header = T)
metadata <- metadata[,-c(12:27)] # remove empty columns

# Change column names
colnames(metadata) <- c("file","file_doi","pub_doi","pub_name","data_doi","first_author_name","year",
                        "species","population","data_type","tissue")
metadata$tissue <- pure_names # switch names back over to Karl's abbreviations

# Load in files
raw_files <- list.files("/workdir/mbb262/kremling_raw_count_v4_hapmap321taxaid/", pattern = ".csv", full.names = TRUE)
raw_files_2 <- list.files("/workdir/mbb262/kremling_raw_count_v4_hapmap321taxaid/", pattern = ".csv")
pure_names <- gsub("_kremling_formatted_v4_hapmapids", "", raw_files_2)
pure_names <- gsub(".csv", "", pure_names)


# Iterate through phenotypes, append tissue name to gene, add metadata,
#   combine to make one large metadata file

# Loop through all 7 tissue files
together <- data.frame()
for (i in seq_len(length(raw_files))){
  # Read in file
  raw_phenos <- data.table::fread(
    input = raw_files[i]
  )

  # Get traits/genes
  raw_phenos <- colnames(raw_phenos)
  raw_phenos <- raw_phenos[-1] # minus out taxa column

  # Paste tissue name in front of gene ID
  raw_phenos <- paste0(pure_names[i], "_", raw_phenos)

  # Add metadata
  meta_gene <- cbind(metadata[i,], data.frame(raw_phenos))

  # Add to outside dataframe
  together <- rbind(together, meta_gene)
}

# Change column name
colnames(together)[12] <- "traits"

# export
write.csv(together, "goodman_kremling_2018_metadata.csv", row.names = F, quote = F)


# -----------------------------------------
# Format trait IDs in fast association data
# (gene ids in this case)
# -----------------------------------------

# TO run on command line
# Rscript /home/mbb262/git_projects/haplotype_v_snp_gwas/src/R/askdb_prep/calc_summary_stats.R > /workdir/mbb262/results_goodman_panel_kremling_formattedtaxaIDs/logging.txt

## Get files ----
assoc_files <- list.files("/workdir/mbb262/results_goodman_panel_kremling/", pattern = ".csv", full.names = TRUE)
assoc_files_2 <- list.files("/workdir/mbb262/results_goodman_panel_kremling/", pattern = ".csv")
assoc_files_2 <- gsub(".csv", "", assoc_files_2)
assoc_files_2

pure_names <- list.files("/workdir/mbb262/results_goodman_panel_kremling/", pattern = ".csv")
pure_names <- gsub("_Kremling_2018_processed.csv", "", pure_names)
pure_names <- gsub(".*_", "", pure_names)
pure_names

## Iterate ----
for (i in seq_len(length(assoc_files))) {
  logging::loginfo(glue::glue("Formatting data for: {assoc_files_2[i]}"))
  logging::loginfo(glue::glue("Using ID for: {pure_names[i]}"))
  
  ## (0) Load in data ----
  ### Association data
  logging::loginfo("Loading in data...")
  assoc_results <- data.table::fread(
    input = assoc_files[i]
  )
  
  # Add on tissue ID
  assoc_results$Trait <- paste0(pure_names[i], "_", assoc_results$Trait)
  
  ## (3) Save to disk ----
  logging::loginfo("Saving to disk...")
  name <- paste0("/workdir/mbb262/results_goodman_panel_kremling_formattedtaxaIDs/", assoc_files_2[i], "_traitsTissueTaxaID.csv")
  data.table::fwrite(
    x = assoc_results, 
    file = name
  )
}
