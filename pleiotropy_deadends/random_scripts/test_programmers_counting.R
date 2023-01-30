# Make example/test data for programmers to test pipeline with


# Load in GWAS results
for (i in 1:10){
  results <- data.table::fread(paste0("~/Desktop/un_subsampled/chrom_", i, "_buckler_2009_processed_00001_filtered.csv")) %>% 
    filter(trait %in% c("DTS_BLUP.Buckler_2009", "ASI_BLUP.Buckler_2009"))
  
  # Subsample
  results <- sample_n(results, 20)
  
  # export
  data.table::fwrite(results, paste0("~/Desktop/subsample test programmers/chrom_", i, "_subsample.csv"))
}

# run script in tassel