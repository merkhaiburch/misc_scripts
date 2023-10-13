# Parse samtools stat function on max reelGene B73 transcripts to all other nam founders

# Packages
library(dplyr)
library(ggplot2)
library(tidyr)

# find all filessas
allFiles <- list.files("/workdir/mbb262/nam_hybrid_rnaseq/references/nam_aligned_transcriptomes/stats", 
                       pattern = ".txt", full.names = T)

# Load all files into R
allStats <- lapply(allFiles, data.table::fread, header = F) %>% data.table::rbindlist(idcol = TRUE)

# Add in identifiers, turn into data frame
shortNames <- gsub("/workdir/mbb262/nam_hybrid_rnaseq/references/nam_aligned_transcriptomes/stats/", "", allFiles)
shortNames <- gsub("_axSplice_eqx_SNsection.txt", "", shortNames)
shortNames <- gsub("_mecatErrorCorrected.contigs", "", shortNames)
shortNames <- data.frame(shortNames)
shortNames$id <- rownames(shortNames) %>% as.numeric()

# Merge together
allStats <- merge(allStats, shortNames, by.x = ".id", by.y = "id") %>% select(-".id")

# Cast data
cast_allStats <- tidyr::pivot_wider(allStats, id_cols = c(shortNames), names_from = V1,
                                    values_from = V2) %>% 
  select("shortNames", "reads mapped", "reads unmapped", "error rate")

# Calculate mapping percentage
cast_allStats$mapping_percentage <- (cast_allStats$`reads mapped`/c(cast_allStats$`reads mapped` + cast_allStats$`reads unmapped`))*100

# Make a ggplot histogram
transcript_mapping_rate <- ggplot(cast_allStats, aes(x = mapping_percentage))+
  geom_histogram(bins = 23)+
  xlab("Mapping Percentage") +
  ylab("Count of NAM Founders")
ggsave("/home/mbb262/git_projects/te_ase_nam/figs/histogram_maxReelGeneB73Transcript_toNAMFounders.png", transcript_mapping_rate)
