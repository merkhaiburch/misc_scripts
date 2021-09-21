# inefficient code to take means within intervals


# Set external variabes
# Set the number of cores to use
numCores <- 3

# Loop through each chromosome, then in parallel analyze the mean LD score
for (j in seq_len(10)){

  # Provide reassuring message on what chromosome we're analyzing
  message(paste0("I am on chromosome: ", j))

  # Subset range file down to just the chromosome we're analyzing
  subset_range <- ranges %>% filter(seqid == j)

  # Load in TASSEL analyzed LD file
  ld_one_chrom <- data.table::fread(input = paste0("/workdir/mbb262/ld_matrices/ld_chr", j, "_282.txt"),
                                    select = c("Locus1", "Position1", "R^2"))

  # Run the parallel statement
  avg_ld_ranges <- c()
  temp <- parallel::mclapply(seq_len(nrow(subset_range)), mc.cores = numCores, function(i){

    logging::loginfo(glue::glue("Mean of: {subset_range$rr_id[i]}"))

    # subset LD file with range intervals
    subset_ld_to_range <- ld_one_chrom[Locus1 == subset_range$seqid[i] & Position1 >= subset_range$start[i] & Position1 <= subset_range$end[i],
                                       .(average_r2_ld = mean(`R^2`))]

    # If LD range returns nothing (i.e. no SNPs averaged over to obtain a avg LD score), provide back a 0
    if (nrow(subset_ld_to_range) == 0){
      subset_ld_to_range <- 0
    }

    # Calculate mean R^2 for this interval
    mean_ld <- data.frame(subset_range[i,], subset_ld_to_range)
    colnames(mean_ld) <- c("seqid", "start", "end", "rr_id", "average_r2_ld")

    # Add results to outside of the loop
    avg_ld_ranges <- rbind(avg_ld_ranges, mean_ld)
    return(avg_ld_ranges)

  })

  # Turn into a dataframe
  temp <- temp %>% rbindlist()

  # Write to file
  data.table::fwrite(temp, file = paste0("/workdir/mbb262/ld_matrices/mean_ld/mean_ld_282_chrom", j, ".csv"))
}


# -----------------------------
# Version 2 but with lapply
# -----------------------------

for (j in seq_len(10)){

  # Provide reassuring message on what chromosome we're analyzing
  message(paste0("I am on chromosome: ", j))

  # Subset range file down to just the chromosome we're analyzing
  subset_range <- ranges %>% filter(seqid == j)

  # Load in TASSEL analyzed LD file
  ld_one_chrom <- data.table::fread(input = paste0("/workdir/mbb262/ld_matrices/ld_chr", j, "_282.txt"),
                                    select = c("Locus1", "Position1", "R^2"))

  # Run the parallel statement
  avg_ld_ranges <- c()
  temp <- lapply(seq_len(nrow(subset_range)), function(i){

    logging::loginfo(glue::glue("Mean of: {subset_range$rr_id[i]}"))

    # subset LD file with range intervals
    subset_ld_to_range <- ld_one_chrom[Locus1 == subset_range$seqid[i] & Position1 >= subset_range$start[i] & Position1 <= subset_range$end[i],
                                       .(average_r2_ld = mean(`R^2`))]

    # If LD range returns nothing (i.e. no SNPs averaged over to obtain a avg LD score), provide back a 0
    if (nrow(subset_ld_to_range) == 0){
      subset_ld_to_range <- 0
    }

    # Calculate mean R^2 for this interval
    mean_ld <- data.frame(subset_range[i,], subset_ld_to_range)
    colnames(mean_ld) <- c("seqid", "start", "end", "rr_id", "average_r2_ld")

    # Add results to outside of the loop
    avg_ld_ranges <- rbind(avg_ld_ranges, mean_ld)
    return(avg_ld_ranges)

  })

  # Turn into a dataframe
  temp <- temp %>% rbindlist()

  # Write to file
  data.table::fwrite(temp, file = paste0("/workdir/mbb262/ld_matrices/mean_ld/mean_ld_282_chrom", j, ".csv"))
}


# do on the tassel test data
# Provide reassuring message on what chromosome we're analyzing
message(paste0("I am on chromosome: ", j))

# Subset range file down to just the chromosome we're analyzing
subset_range <- ranges %>% filter(seqid == j)

# Load in TASSEL analyzed LD file
ld_one_chrom <- data.table::fread(input = "~/Desktop/ld_tassel_output.txt",select = c("Locus1", "Position1", "R^2"))

# Run the parallel statement
avg_ld_ranges <- c()
temp <- lapply(seq_len(nrow(ranges)), function(i){

  logging::loginfo(glue::glue("Mean of: {ranges$rr_id[i]}"))

  # subset LD file with range intervals
  subset_ld_to_range <- ld_one_chrom[Locus1 == ranges$seqid[i] & Position1 >= ranges$start[i] & Position1 <= ranges$end[i],
                                     .(average_r2_ld = mean(`R^2`))]

  # If LD range returns nothing (i.e. no SNPs averaged over to obtain a avg LD score), provide back a 0
  if (nrow(subset_ld_to_range) == 0){
    subset_ld_to_range <- 0
  }

  # Calculate mean R^2 for this interval
  mean_ld <- data.frame(ranges[i,], subset_ld_to_range)
  colnames(mean_ld) <- c("seqid", "start", "end", "rr_id", "average_r2_ld")

  # Add results to outside of the loop
  avg_ld_ranges <- rbind(avg_ld_ranges, mean_ld)
  return(avg_ld_ranges)

})

# Turn into a dataframe
temp <- temp %>% rbindlist()

# Write to file
data.table::fwrite(temp, file = paste0("/workdir/mbb262/ld_matrices/mean_ld/mean_ld_282_chrom", j, ".csv"))

data.table::fwrite(temp, file = "~/Downloads/example_ld_mean_output.txt")
