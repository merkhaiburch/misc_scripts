# ---------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-04-12 
# Updated... 2022-04-12

# Description 
# Count the number of GWA hits in two interval types:
# - all B73 v5 transcripts from start to stop
# - all B73 v5 transcripts from start to stop +5kb window up and downstream
# --------------------------------------------------------------------------

# Load helpful packages -----------------------------------------------------

library(dplyr)
library(data.table)
library(GenomicRanges)


# Gather GWAS results from blfs1 --------------------------------------------

# field NAM
# scp -r mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/mbb262/results/pleiotropy/gwa_rtassel_all_nam_traits/nam_all/processed_results/* /workdir/mbb262/results/nam


# Create intervals -----------------------------------------------------------

# Intervals from:
# cd /workdir/mbb262
# https://datacommons.cyverse.org/browse/iplant/home/maizegdb/maizegdb/B73v5_JBROWSE_AND_ANALYSES/B73v3-B73v5_and_B73v4-B73v5_gene_model_associations/B73v4_B73v5_liftoff_genemodel_CDS_xref.txt.zip

# Load in interval file
intervals <- read.table("~/Downloads/B73v4_B73v5_liftoff_genemodel_CDS_xref.txt") %>% 
  select(-"V5", -"V11")

# Change column names
colnames(intervals) <- c("v4_chrom", "v4_start", "v4_end", "v4_gene", "v4_strand", 
                         "v5_chrom", "v5_start", "v5_end", "v5_gene", "v5_strand",
                         "v5_length")
intervals$v4_chrom <- gsub("chr", "", intervals$v4_chrom)
intervals$v4_chrom <- as.numeric(as.character(intervals$v4_chrom))

# Remove v4 genes that do not map to v5, and any v5 genes that don't map to v4
intervals <- intervals %>% filter(v5_gene != ".")
intervals <- intervals[complete.cases(intervals), ]

# make file with chrom, start, stop --> no window around genes
intervals_1 <- intervals %>% select("v4_chrom", "v4_start", "v4_end", "v4_gene", "v4_strand")

# make file with chrom, start, stop --> 5 kb window around genes
intervals_2 <- GenomicRanges::makeGRangesFromDataFrame(intervals_1, strand.field = "v4_strand")
temp <- resize(intervals_2, width(intervals_2)-10, fix = "both")
# temp <- resize(temp, width(temp)+10, fix = "end")

intervals_2 # 34722-38366
temp # 34712-38376

intervals_2 <- intervals %>% select("v4_chrom", "v4_start", "v4_end", "v4_gene")
temp <- intervals$v4_start-intervals$v4_end

# Export without windows around genes
write.csv(intervals_1, "/workdir/mbb262/aimee_intervals_1.csv", row.names = FALSE, quote = FALSE)
write.csv(intervals_2, "/workdir/mbb262/aimee_intervals_2.csv", row.names = FALSE, quote = FALSE)


# Run Zack's counting pipeline ------------------------------------------------

# find all the files we need
system("find /workdir/mbb262/results/nam/*.csv -type f > /workdir/mbb262/nam_files.list")
system("find /workdir/mbb262/results/goodman/*.csv -type f > /workdir/mbb262/goodman_files.list")

# Make directories for results --> interval 1 (no windows)
system("mkdir -k /workdir/mbb262/output_counts_interval_1/nam")
system("mkdir /workdir/mbb262/output_counts_interval_1/goodman")

# Make directories for results --> interval 2 (5kb windows)
system("mkdir -k /workdir/mbb262/output_counts_interval_2/nam")
system("mkdir /workdir/mbb262/output_counts_interval_2/goodman")

# run count association plugin in parallel --> interval_1 files (no window around genes)
system("/programs/parallel/bin/parallel -j 20 '/home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl -debug /workdir/mbb262/output_counts_interval_1/nam/debug.txt -Xmx200g -CountAssociationsPlugin -intervals /workdir/mbb262/aimee_intervals_1.csv -gwasResults {} -outputCountFile /workdir/mbb262/output_counts_interval_1/nam/outputCount_{/.}.txt -endPlugin' :::: /workdir/mbb262/results/nam/nam_files.list")
system("/programs/parallel/bin/parallel -j 20 '/home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl -debug /workdir/mbb262/output_counts_interval_1/goodman/debug.txt -Xmx200g -CountAssociationsPlugin -intervals /workdir/mbb262/aimee_intervals_1.csv -gwasResults {} -outputCountFile /workdir/mbb262/output_counts_interval_1/goodman/outputCount_{/.}.txt -endPlugin' :::: /workdir/mbb262/results/goodman/goodman_files.list")

# run count association plugin in parallel --> interval_1 files (no window around genes)
system("/programs/parallel/bin/parallel -j 20 '/home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl -debug /workdir/mbb262/output_counts_interval_2/nam/debug.txt -Xmx200g -CountAssociationsPlugin -intervals /workdir/mbb262/aimee_intervals_2.csv -gwasResults {} -outputCountFile /workdir/mbb262/output_counts_interval_2/nam/outputCount_{/.}.txt -endPlugin' :::: /workdir/mbb262/results/nam/nam_files.list")
system("/programs/parallel/bin/parallel -j 20 '/home/mbb262/bioinformatics/tassel-5-standalone/run_pipeline.pl -debug /workdir/mbb262/output_counts_interval_2/goodman/debug.txt -Xmx200g -CountAssociationsPlugin -intervals /workdir/mbb262/aimee_intervals_2.csv -gwasResults {} -outputCountFile /workdir/mbb262/output_counts_interval_2/goodman/outputCount_{/.}.txt -endPlugin' :::: /workdir/mbb262/results/goodman/goodman_files.list")


# Function to aggregate Zack's counts by chromosome ----------------------------

# Zack's count function in TASSEL returns a table with all intervals and their 
# counts of GWA hits, my inputs were chrom by chrom gwa hits meaning that 1 
# output file would contain count results for 1 chromosome and 0 counts for 
# the other 9 chromosomes, this function drops those empty 9 chromosomes, 
# combines all 10 chroms worth of results, and then exports that dataframe
aggregate_gwa_counts <- function(data_path, out_dir, data_study, pattern_file_end, name_append) {
  # data_path = directory where results from Zack's count function reside
  # out_dir = directory where you want aggregated results to live after analyzing
  # data_study = a vector containing the order of the study name of each file from the data_path directory
  
  # Unique studies to iterate by
  unique_study <- unique(data_study)
  
  # Keep track of all data in this path's directory
  data_all_files <- list.files(data_path, pattern = pattern_file_end)
  
  # Use lapply loop
  for (i in seq_along(unique_study)){
    
    # Make an empty list to aggregate results
    temp <- data.table()
    
    # status message
    message("On trait: ", unique_study[i])
    
    # Iterate through all 10 chromosomes worth of a single study's results
    for (j in seq_len(10)){
      
      # Load in interval file, subset to sepcific chromosome
      intervals_b73 <- data.table::fread("/workdir/mbb262/genic_intergenic_intervals_b73_v4.49.csv") %>% 
        filter(seqnames == j) %>% 
        select(seqnames, rr_id)
      
      # Extract files matching the name of a single study 
      full_names <- data_all_files[grep(unique_study[i], x = data_all_files, fixed = TRUE)]
      
      # Just select the file with the matching j chromosome
      full_names <- full_names[grep(paste0("chrom_", j, "_"), full_names)]
      # print(full_names)
      
      # Read in the file with fread
      results <- data.table::fread(paste0(data_path, full_names), nThread = n_threads)
      
      # only keep rows (i.e. intervals on this specific chromosome)
      subset_results <- results[intervals_b73, on = .(Interval = rr_id)] 
      
      # Rearrange seqnames column
      subset_results <- subset_results[,seqnames:=NULL]
      
      # combine with results outside of this lapply statement
      temp <- bind_rows(temp, subset_results)
      
    }
    # Save to file
    data.table::fwrite(temp, file = paste0(out_dir, unique_study[i], name_append), nThread = n_threads)
    print(dim(temp))
  }
}

# Function to aggregate all trait data sets together
aggregate_trait <- function(file_path, output_id) {
  # Load in files under the specified path
  file_files <- list.files(path = file_path, full.names = TRUE)
  
  # Load in everything with a match
  file_to_agg <- Reduce(function(x, y) merge(x, y, by = "Interval"), 
                        lapply(file_files, data.table::fread, nThread = n_threads))
  
  # Write to file
  data.table::fwrite(file_to_agg, file = output_id, nThread = n_threads)
  print(dim(file_to_agg))
}


# Use functions on different datasets ------------------------------------------

# Make result directories
system("mkdir -p /workdir/mbb262/aggregated_counts/nam_filtered_count_agg/")

# NAM --> interval 1
nam_interval_1_path <- "/workdir/mbb262/output_counts_interval_1/nam/"
nam_interval_1_study <- list.files(nam_interval_1_path, pattern = "*_outputCount.txt", include.dirs = FALSE)
nam_interval_1_study <- gsub("_peaks.csv_outputCount.txt", "", nam_interval_1_study)
nam_interval_1_study <- gsub("chrom_[0-9]+_", "", nam_interval_1_study)

aggregate_gwa_counts(data_path = nam_interval_1_path,
                     out_dir = "/workdir/mbb262/aggregated_counts/nam_filtered_count_agg/",
                     data_study = nam_interval_1_study,
                     pattern_file_end = "*_outputCount.txt",
                     name_append = "aggregatedCount10Chroms_nam.txt")


# Nam filtered all --> DONE 2022-01-29
nfa_path <- "/workdir/mbb262/aggregated_counts/nam_filtered_count_agg/"
aggregate_trait(file_path = nfa_path, output_id = "nam_filtered_all_traits.txt")


# Count number of traits across rows --------------------------------------------

# Datatbale function to replace NAs with 0's
replace_nas <- function(DT) {
  for (j in seq_len(ncol(DT)))
    set(DT,which(is.na(DT[[j]])),j,0)
  
  return(DT)
}

# NOTE: script 18_aggregate_gwa_counts leaves NAs when there are no GWA hits in any intervals
# we need to replace those NAs with 0's (e.g. there are only 3 sig. snps for a trait, the focal
# interval will have a count of 3 but all other intervals will have NA values)

# Start loading in data
data_path <- "/workdir/mbb262/aggregated_counts/aggregated_by_pop_or_trait/"

# - Nam filtered physiological (nfp)
nfp <- data.table::fread(paste0(data_path, "nam_filtered_all_traits.txt"))
nfp[is.na(nfp)] <- 0

# Nam filtered physiological (nfp)
nfp_pleiotropy <- calc_pleiotropy_countUnique(nfp, "nam_filtered_all")
