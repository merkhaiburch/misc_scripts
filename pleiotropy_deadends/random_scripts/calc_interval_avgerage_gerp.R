# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch, Emily Yi
# Contact... mbb262@cornell.edu, ety8@cornell.edu
# Date...... 2021-07-16
#
# Description 
#   - average GERP scores in intervals
# ---------------------------------------------------------------


# ----------------------------
# Load in helpful packages
# -----------------------------
library(dplyr)
library(data.table)
library(parallel)
library(tidyr)
library(assertthat)


# -----------------------------
# Gather genic/intergenic intervals
# -----------------------------
# Get intervals from cbsu
# scp mbb262@cbsublfs1.tc.cornell.edu:/data1/users/mbb262/haplotype_ranges/genic_intergenic_intervals_b73_v4.49.csv /workdir/mbb262
### Genic and intergenic ranges (not from the phg but from genic and intergenic ranges in v4 b73)
# these are merged intervals

ranges <- data.table::fread("/home/ety8/gerp/genic_intergenic_intervals_b73_v4.49.csv")

# Change column names
colnames(ranges) <- c("seqid", "start", "end", "rr_id")

# New column for average GERP
ranges$avg_gerp <- rep(20, nrow(ranges))


# -------------------
# Gather GERP data
# -------------------
# gerp <- data.table::fread("/workdir/mbb262/Zea_mays.allChr.rates", header = FALSE)
# fixed spacing of above file with sed, saved as gerp_tab.txt
gerp <- fread("/home/ety8/gerp/gerp_tab.txt")
# gerp <- fread("/home/ety8/gerp/gerp_sample.txt") # this was used for quick debugging
assert_that(all(complete.cases(gerp))) # checks for NA, NaN
colnames(gerp) <- c("chr", "pos", "neutral_rate", "gerp_score")


# ----------------------
# Start analyzing data
# ----------------------

# Set external variabes
# Set the number of cores to use
numCores <- 3

# Loop through each chromosome, then in parallel analyze the mean GERP score
for (chrom in seq_len(10)){
  
  # Provide reassuring message on what chromosome we're analyzing
  message(paste0("I am on chromosome: ", chrom))
  
  # Subset range file down to just the chromosome we're analyzing
  subset_range <- ranges[seqid == chrom]
  subset_gerp_chrom <- gerp[chr == chrom]
  
  # Run the parallel statement
  for (i in seq_len(nrow(subset_range))) {

    logging::loginfo(glue::glue("Mean of: {subset_range$rr_id[i]}"))
    
    # subset GERP file to range interval
    subset_gerp_to_range <- subset_gerp_chrom[pos >= subset_range$start[i] & pos <= subset_range$end[i],
                                 gerp_score]
    
    # Calculate mean GERP for this interval
    if (subset_gerp_to_range %>% length == 0)
      mean_gerp <- NA_real_
    else {
      mean_gerp <- mean(subset_gerp_to_range)
      assert_that(!is.nan(mean_gerp))
    }
    
    # return values
    ranges[seqid == chrom & rr_id == subset_range$rr_id[i],
           avg_gerp := mean_gerp]
  }
}


# Write to file
data.table::fwrite(ranges, file = "avg_gerp.csv")


saved <- fread("/home/ety8/gerp/avg_gerp.csv")
assert_that(!any(na.omit(saved$avg_gerp) == 20)) # check that there are no 20s left

# Ideally, all genic regions have gerp (not NA or NaN) - but maybe acceptable