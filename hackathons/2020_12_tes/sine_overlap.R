# ---------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2020-12-16 
#
# Description 
#   - December 2020 hackathon: TE annotation group
#   - Minimap known SINE TEs to maize v5, compare with sine-finder results
# ---------------------------------------------------------------

# Load packages
library(GenomicRanges)
library(Rsamtools)
library(data.table)
library(dplyr)

# Load in sam file with found SINEs in V5
bam <- Rsamtools::scanBam("~/Downloads/sine_alignment.bam")

# Turn into a dataframe
bam <- do.call("DataFrame", bam) %>% as.data.table()
bam$end <- bam$pos+bam$qwidth
bam$where <- rep("minimap", nrow(bam))

# turn into genomic range object
bam <- GenomicRanges::makeGRangesFromDataFrame(bam, keep.extra.columns = TRUE, 
                                                start.field = "pos", end.field = "end",
                                                seqnames.field = "rname")

# Load in SINE finder results in gff file
sf <- ape::read.gff("~/Downloads/B73v5.RST.stringentAB.gff3")
sf$where <- rep("sineFinder", nrow(sf))
sf <- GenomicRanges::makeGRangesFromDataFrame(sf, keep.extra.columns = TRUE)

# Use genomicRanges to find overlaps between the two files
temp <- subsetByOverlaps(bam, sf) %>% as.data.frame()
temp$overlap <- width(intersect(bam, sf))
table(temp$qname) # count # overlaps by family

# What overlapped in sf?
temp2 <- subsetByOverlaps(sf, bam) %>% as.data.frame()
bam_sf_overlap <- GenomicRanges::makeGRangesFromDataFrame(temp2, keep.extra.columns = TRUE)

# What fspecific familes were found with sinefinder
tmp <- read.table("~/Downloads/quick_sine_stringentAB_famclust.txt", header = TRUE)
tmp$chrom <- gsub("_.*", "", tmp$loc)
tmp$pos <- gsub("chr[0-9]{1,2}_R_", "", tmp$loc)
tmp$pos <- gsub("chr[0-9]{1,2}_S_", "", tmp$pos)
tmp$pos <- gsub("chr[0-9]{1,2}_F_", "", tmp$pos)
tmp$pos <- gsub("_.*", "", tmp$pos)
tmp <- tmp %>% select(-loc) %>% filter(chrom != "RST")
tmp$start <- gsub(":.*", "", tmp$pos) %>% as.numeric(as.character())
tmp$end <- gsub("*.:", "", tmp$pos) %>% as.numeric(as.character())
tmp <- GenomicRanges::makeGRangesFromDataFrame(sf, keep.extra.columns = TRUE)

# of the SINES that overlap with minimap2 sines, how many were ID'd with clustering?
subsetByOverlaps(tmp, bam_sf_overlap) %>% as.data.frame()
countOverlaps(tmp, bam_sf_overlap)


# V1 sine bam
# Load in sam file with found SINEs in V5
bam_v1 <- Rsamtools::scanBam("~/Downloads/v1_sine_alignment.bam")

# Turn into a dataframe
bam_v1 <- do.call("DataFrame", bam_v1) %>% as.data.table()
bam_v1$end <- bam_v1$pos+bam_v1$qwidth
bam_v1$where <- rep("minimap", nrow(bam_v1))

# Make table of counts by family
table(bam_v1$qname)

# turn into genomic range object
bam_v1 <- GenomicRanges::makeGRangesFromDataFrame(bam_v1, keep.extra.columns = TRUE, 
                                               start.field = "pos", end.field = "end",
                                               seqnames.field = "rname")


# Ds1 alignmenys

# Ds1 against v5
bam_v5_ds1 <- Rsamtools::scanBam("~/Downloads/ds1_alignment.bam")
bam_v5_ds1 <- do.call("DataFrame", bam_v5_ds1) %>% as.data.table()
dim(bam_v5_ds1)

bam_v2_ds1 <- Rsamtools::scanBam("~/Downloads/v2_ds1_alignment.bam")
bam_v2_ds1 <- do.call("DataFrame", bam_v2_ds1) %>% as.data.table()
dim(bam_v2_ds1)

bam_v1_ds1 <- Rsamtools::scanBam("~/Downloads/v1_ds1_alignment.bam")
bam_v1_ds1 <- do.call("DataFrame", bam_v1_ds1) %>% as.data.table()
dim(bam_v1_ds1)
