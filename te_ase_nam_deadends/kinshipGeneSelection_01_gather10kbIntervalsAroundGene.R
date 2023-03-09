# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-01-09
# Updated... 2023-01-24
#
# Description:
# Get a 100 kb interval around a gene, make bed file
# ------------------------------------------------------------------------------

# Libraries
library(magrittr)
library(dplyr)
library(GenomicRanges)

## Test that this windowing method works in genomic ranges ---------------------
updownwindow <- 10
temp <- data.frame(A = c(10,20,30), B = c(100,200,300), seqid = c(1,2,3), strand = c("+", "-", "+"))
gr_temp <- GenomicRanges::makeGRangesFromDataFrame(df = temp, 
                                                   start.field = "A", 
                                                   end.field = "B", 
                                                   seqnames.field = "seqid")
gr_temp <- resize(gr_temp, width=width(gr_temp) + updownwindow, fix="start")
gr_temp <- resize(gr_temp, width=width(gr_temp) + updownwindow, fix="end")
gr_temp

# Merge any overlapping ranges
reduce(gr_temp)

# Turn back into R object
lala <- data.frame(gr_temp)
lala


## Get intervals for v5 genes --------------------------------------------------
# Get the start and stop positions 10 kb up and downstream of each gene

# Load in v5 gff file
# https://www.maizegdb.org/genome/assembly/Zm-B73-REFERENCE-GRAMENE-4.0#downloads
# Get all maize v4 genes from gff file 
gff <- ape::read.gff("~/Downloads/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3", 
                     na.strings = c(".", "?"), GFF3 = TRUE) %>% 
  filter(type == "gene")

# Remove extra columns and scaffolds
gff <- gff %>% 
  select("seqid", "start", "end", "strand") %>% 
  filter(seqid %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10"))

# Turn into a genomic ranges object
gr_gff <- GenomicRanges::makeGRangesFromDataFrame(df = gff)

# Add on 10 kb windows
window_size <- 10000
gr_gff <- resize(gr_gff, width=width(gr_gff) + window_size, fix="start")
gr_gff <- resize(gr_gff, width=width(gr_gff) + window_size, fix="end")

# Merge any overlapping ranges
gr_gff <- reduce(gr_gff)
temp <- width(gr_gff)
summary(temp)

# Turn back into df
df_gr_gff <- data.frame(gr_gff)

# If the start position is negative, replace with 1 (becomes 0 when turned into bed object)
df_gr_gff$start[df_gr_gff$start < 0] <- 1

# Turn back to genomic ranges object
gr <- GenomicRanges::makeGRangesFromDataFrame(df = df_gr_gff)

# Export bed file
df <- data.frame(seqnames=seqnames(gr),
                 starts=start(gr)-1,
                 ends=end(gr))

write.table(df, file="~/git_projects/te_ase_nam/data/10kb_up_down_gene.bed", 
            quote=F, sep="\t", row.names=F, col.names=F)




