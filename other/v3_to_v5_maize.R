# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2021-06-21 
#
# Description 
#   - converting from v3 to v5 gene ids
#   - 
# ---------------------------------------------------------------

# Get data from cyverse
# https://datacommons.cyverse.org/browse/iplant/home/maizegdb/maizegdb/B73v5_JBROWSE_AND_ANALYSES/B73v3-B73v5_and_B73v4-B73v5_gene_model_associations/B73v3_B73v5_liftoff_genemodel_CDS_xref.txt.zip

# Get Anju's translation file
# Slack

# Load packages
library(magrittr)
library(dplyr)


# -----------------
# Load data into R
# -----------------

# cross reference file
ref <- read.table("B73v3_B73v5_liftover_genemodel_CDS_xref.txt", 
                  sep = "\t", na.strings = c("."), header = F)

colnames(ref) = c("v3_chrom", "v3_start", "v3_end", "v3_geneid", "v3_gene_id2", "v3_strand",
                  "v5_chrom", "v5_start", "v5_end", "v5_gene_id", "v5_gene_id2", "v5_strand",
                  "v5_width")

# Remove extra columns
ref <- ref |> select("v3_chrom", "v3_start", "v3_end", "v3_geneid", 
                     "v5_chrom", "v5_start", "v5_end", "v5_gene_id")

# anju file
hare <- read.table("genelist_included_HARE.csv", header = TRUE, sep = ",", row.names = 1)
colnames(hare) <- "hare_v3"

# Merge cross ref table with hare IDs
tov5 <- merge(x = ref, y = hare,
              by.x = "v3_geneid", by.y = "hare_v3")

# Format chromosome names
tov5$v5_chrom <- gsub("chr", "", tov5$v5_chrom) |> as.numeric(as.character())

# Only save unique gene ids from v3
tov5 <- distinct(tov5, v3_geneid, .keep_all = TRUE)

# keep only columns we need for bed file
tov5 <- tov5 |> select("v5_chrom", "v5_start", "v5_end") |> 
  arrange(v5_chrom, v5_start)

# Remove rows with NA values
tov5 <- tov5[complete.cases(tov5), ]

# Save as a bed file
write.table(tov5, file = "hare_v5_intervals.bed", 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)



