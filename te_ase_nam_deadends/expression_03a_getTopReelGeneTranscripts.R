# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-06-27
# Updated... 2023-06-27
#
# Description:
# Select the best gene model based off reelGene scores
# https://bitbucket.org/bucklerlab/reelgene_scores/src/master/
# ------------------------------------------------------------------------------

# Load packages
library(dplyr)
library(tidyr)

# Load in table with scores
reel <- read.delim("~/git_projects/reelgene_scores/reelGeneScores_release_20230314.txt")

# Add a column with just the base gene name and transcript identifier
reel <- reel %>%
  separate(transcript, c("gene", "transcriptID"), "_")

# Sum up all four scores
reel$sumCols <- rowSums(reel[,c("reelProteinScore", "reelExon", "trscStartJunctionScore", "trscStopJunctionScore")], na.rm=TRUE)

# Select the gene's transcript with the highest score 
top_reel <- reel %>% 
  group_by(gene) %>% 
  slice(which.max(sumCols))

# Remake into gene_transcriptID
top_genes <- data.frame(paste0(top_reel$gene, "_", top_reel$transcriptID))
colnames(top_genes) <- "transcriptID"

# Export to file
write.csv(top_genes, 
          file = "~/git_projects/te_ase_nam/data/top_reel_transcripts.txt", 
          quote = FALSE, 
          row.names = FALSE)
