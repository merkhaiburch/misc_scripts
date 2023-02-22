# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-01-25
# Updated... 2023-01-25
#
# Description:
# Create an interval file that contains the interval up until the next gene
# upstream and downstream. 
# ------------------------------------------------------------------------------

# Libraries
library(magrittr)
library(dplyr)
library(GenomicRanges)


# Proof of concept that my neighboring interval method works -------------------

# Create dummy data
# Threw in two weird cases, one that overlaps another range, one that's completely within another
temp <- data.frame(seqnames = c(1,1,1,1,1,1), 
                   strand = c("+", "+", "+", "+", "+", "-"),
                   start = c(3,9,20,35,40,40), 
                   end = c(10,15,30,50,60,50), 
                   ID = c("A","B", "C", "D", "E", "F"))
temp

# Turn into GRanges object
temp_gr <- GenomicRanges::makeGRangesFromDataFrame(df = temp, keep.extra.columns = TRUE)
temp_gr

# Get intervals in gaps
temp_gr_intergenic <- GenomicRanges::gaps(temp_gr)
temp_gr_intergenic
temp_gr_intergenic_df <- data.frame(temp_gr_intergenic) %>% select(-width, -strand, -seqnames)

# Find flanking intervals for each gene
temp_intervals_2_right <- GenomicRanges::precede(temp_gr, temp_gr_intergenic, select = "all") %>% data.frame()
temp_intervals_2_right$direction <- rep("right", nrow(temp_intervals_2_right))

temp_intervals_2_left <- GenomicRanges::follow(temp_gr, temp_gr_intergenic) %>% data.frame()
temp_intervals_2_left$queryHits <- rownames(temp_intervals_2_left)
colnames(temp_intervals_2_left)[1] <- "subjectHits"
temp_intervals_2_left$direction <- rep("left", nrow(temp_intervals_2_left))

temp_flank <- rbind(temp_intervals_2_right, temp_intervals_2_left[,c(2,1,3)]) %>% arrange(queryHits, subjectHits)
temp_flank

# Bring back in gene names, intervals for subject hits, rename columns
temp_ids <- merge(temp, temp_flank, by.x = 0, by.y = "queryHits", all = TRUE) %>% 
  select(-"Row.names")
temp_ids <- merge(temp_ids, temp_gr_intergenic_df, by.x = "subjectHits", by.y = 0, all = TRUE) %>% 
  select(-subjectHits)
colnames(temp_ids) <- c("seqnames", "strand", "start_gene", "end_gene", "gene", 
                        "direction", "start_flank", "end_flank")
temp_ids %>% 
  arrange(gene, start_flank)

# Remove variables
rm(temp)
rm(temp_gr)
rm(temp_gr_intergenic)
rm(temp_gr_intergenic_df)
rm(temp_intervals_2_right)
rm(temp_intervals_2_left)
rm(temp_flank)
rm(temp_ids)


## Get intervals for v5 genes --------------------------------------------------

# Load in v5 gff file
# https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz
# Get all maize v4 genes from gff file 
gff <- ape::read.gff("~/Downloads/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3", 
                     na.strings = c(".", "?"), GFF3 = TRUE) 

# Gsub out gene name
gff$gene <- gsub("ID=", "", gff$attributes)
gff$gene <- gsub(";.*", "", gff$gene)

# Get chromosome information
chroms <- gff %>% 
  filter(type == "chromosome", seqid %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10")) %>% 
  select("seqid", "start", "end")

# Remove extra columns and scaffolds
gff <- gff %>% 
  select("seqid", "start", "end", "strand", "gene", "type") %>% 
  filter(seqid %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10"), type == "gene") %>% 
  select(-"type")
rownames(gff) = seq(length=nrow(gff))

# Turn into a gr object
gff_gr <- GenomicRanges::makeGRangesFromDataFrame(df = gff, keep.extra.columns = TRUE)

# Get intervals in gaps
gff_gr_intergenic <- GenomicRanges::gaps(gff_gr)
gff_gr_intergenic
gff_gr_intergenic_df <- data.frame(gff_gr_intergenic) %>% select(-width, -strand, -seqnames)

# Find flanking intervals for each gene
gff_intervals_2_right <- GenomicRanges::precede(gff_gr, gff_gr_intergenic, select = "all") %>% data.frame()
gff_intervals_2_right$direction <- rep("right", nrow(gff_intervals_2_right))

gff_intervals_2_left <- GenomicRanges::follow(gff_gr, gff_gr_intergenic) %>% data.frame()
gff_intervals_2_left$direction <- rep("left", nrow(gff_intervals_2_left))
gff_intervals_2_left$queryHits <- rownames(gff_intervals_2_left)

colnames(gff_intervals_2_left)[1] <- "subjectHits"
gff_flank <- rbind(gff_intervals_2_right, gff_intervals_2_left[,c(2,1,3)]) %>% arrange(queryHits, subjectHits)
head(gff_flank)
summary(gff_flank)

# Bring back in gene names, intervals for subject hits, rename columns
gff_ids <- merge(gff, gff_flank, by.x = 0, by.y = "queryHits", all = TRUE) %>% 
  select(-"Row.names")
gff_ids <- merge(gff_ids, gff_gr_intergenic_df, by.x = "subjectHits", by.y = 0, all = TRUE) %>% 
  select(-subjectHits)
colnames(gff_ids) <- c("seqnames", "start_gene", "end_gene", "strand", "gene", 
                       "direction", "start_flank", "end_flank")
gff_ids <- gff_ids %>% arrange(gene, start_flank)
head(gff_ids)
dim(gff_ids)

# Find rows with NAs, these are the last genes on each chromosome, all are on the negative
# strand, need to create the interval to the left to the end of the chromosome
# Fix the issue and add back to main df
gff_ids[rowSums(is.na(gff_ids)) > 0,] 
ends_gff_ids <- gff_ids[rowSums(is.na(gff_ids)) > 0,] # subset them out from df
ends_gff_ids$start_flank <- ends_gff_ids$end_gene + 1 # start_flank position starts +1 from end (verified this with unit test example)
ends_gff_ids <- merge(ends_gff_ids, chroms, by.x = "seqnames", by.y = "seqid") %>% select(-"end_flank",-"start") # end = end of chromosome
colnames(ends_gff_ids)[8] <- "end_flank"
gff_ids <- gff_ids[complete.cases(gff_ids), ] # remove old lefts from df
gff_ids <- rbind(gff_ids, ends_gff_ids) # add new lefts to the df
dim(gff_ids)

# Find the genes with no right interval. These are all on the positive strand and are
# the last genes on the chromosome arms. Make them a right interval.
temp <- data.frame(table(gff_ids$gene)) %>% filter(Freq == 1) %>% select(Var1)
(ends_pos_gff_ids <- gff_ids%>% filter(gene %in% temp$Var1))
ends_pos_gff_ids <- ends_pos_gff_ids %>% select(-"direction", -"start_flank", -"end_flank")
ends_pos_gff_ids$direction <- rep("right", nrow(ends_pos_gff_ids)) # add "right" direction
ends_pos_gff_ids$start_flank <- ends_pos_gff_ids$end_gene + 1 # start flank is gene's end + 1
ends_pos_gff_ids <- merge(ends_pos_gff_ids, chroms, by.x = "seqnames", by.y = "seqid") %>% select(-"start") # end = end of chromosome
colnames(ends_pos_gff_ids)[8] <- "end_flank"
gff_ids <- rbind(gff_ids, ends_pos_gff_ids) # add new rights to the df


# Check results ----------------------------------------------------------------

# Should have the number of genes *2 number of intervals
(nrow(gff)*2)
nrow(gff_ids)
(nrow(gff)*2)==nrow(gff_ids)

# Table to check that each gene has an interval to the left and right (should == 0)
data.frame(table(gff_ids$gene)) %>% filter(Freq == 1)

# See if any genes have two rights or two lefts instead of one each (should have nothing)
gff_ids %>% group_by(gene, direction) %>% summarise(n = n()) %>% filter(n == 2)

# Get stats
gff_ids$width_intervals <- gff_ids$end_flank - gff_ids$start_flank
summary(gff_ids$width_intervals) # mean = 104909, median = 55727, max = 4858221

# Plot
library(ggplot2)
a <- ggplot(gff_ids, aes(x=width_intervals)) +
  geom_histogram(bins = 100) +
  theme_bw()
ggsave(a, filename = "~/git_projects/te_ase_nam/images/histogram_size_of_intervals.png")

# Export results to file
write.csv(gff_ids, file = "~/git_projects/te_ase_nam/data/v5_gene_flanking_until_next_gene.csv",
          row.names = F, quote = F)



