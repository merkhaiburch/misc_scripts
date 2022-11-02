# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-07-27 
# Updated... 2022-07-27

# Description 
# Calculate and correlate different genome-wide kinship matrices 
# containing maize functional components for use in mixed models
#
# Kinship matrix inputs:
# 1 - K genome wide
# 2 - K 5' UTR
# 3 - K TF binding sites
# 4 - K TF binding sites + UTR
# 5 - K TEs within 2 kb of a gene?
# ---------------------------------------------------------------

# Clear environment
rm(list = ls())
gc()

# Set memory
options(java.parameters = "-Xmx100g")
options(java.parameters = "-Xmx100g")

library(rJava)
memUsageAfter <- .jcall(.jnew("java/lang/Runtime"), "J", "maxMemory") / 1024^3
message("Your max memory usage after changing it is: ", memUsageAfter, "GB")

# Load in packages
library(rTASSEL)
library(dplyr)
library(data.table)

# Start logger
rTASSEL::startLogger(fullPath = "/workdir/mbb262", fileName = "rtassel_logger.txt")


# Work with vcf file -----------------------------------------------------------

# scp -r mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/mbb262/te_ase_nam/genotypes_for_k/v5_inbred_hybrid_genotypes /workdir/mbb262/genotypes

# Load in vcf file
tasGeno <-  rTASSEL::readGenotypeTableFromPath("/workdir/mbb262/genotypes/v5_inbred_hybrid_genotypes/hmp321_282_agpv5_inbred_hybrid_all_chroms.vcf.gz",
                                               keepDepth = FALSE)
print(tasGeno)

# of SNPs: 10156186
# of taxa: 113


# 1 - K genome wide ---------------------------------------------

# Leave vcf file intact
# of SNPs: 10156186

# calculate kinship
tasKin1 <- kinshipMatrix(tasObj = tasGeno)


# 2 - K 5' UTR + fixed window size ---------------------------------------------

# get gff file
system("wget -P /workdir/mbb262/ https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz")

# Load in gff file, filter to 5' UTR regions of canonical transcripts
gff_5utr <- ape::read.gff("/workdir/mbb262/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz", 
                          na.strings = c(".", "?"), GFF3 = TRUE) %>% 
  filter(type == "five_prime_UTR", seqid %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10")) %>% 
  select("seqid", "start", "end", "strand", "attributes")

# Parse out gene names from metadata column & Subset gff just canonical transcripts
system("wget -P /workdir/mbb262/ https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.canonical_transcripts")
gff_5utr$attributes <- gsub("Parent=", "", gff_5utr$attributes)
canonical <- read.delim("/workdir/mbb262/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.canonical_transcripts", header = FALSE)
dim(gff_5utr)
gff_5utr <- gff_5utr %>% filter(attributes %in% canonical$V1) 
dim(gff_5utr)

# add a fixed window size to each interval
window_size <- 500
gff_5utr$start <- gff_5utr$end - window_size

# Turn into a genomic ranges object
gr_5utr <- GenomicRanges::makeGRangesFromDataFrame(gff_5utr, keep.extra.columns=TRUE)

# Subset vcf file using Genomic Ranges object
tasGeno_5utr <- tasGeno %>% 
  rTASSEL::filterGenotypeTableSites(
    siteRangeFilterType = "none",
    gRangesObj = gr_5utr)
# of sites left: 75399

# calculate kinship
tasKin2 <- kinshipMatrix(tasObj = tasGeno_5utr)


# 3 - K TF binding sites -------------------------------------------------------

# File from Tu 2020 nature communications paper, supp. table 5
# https://www.nature.com/articles/s41467-020-18832-8
# File was in v4 coordinates, used Crossmap to go to v5 (see uplift_hapmap_v5.bash script)

# Load in bed file 
tf <- read.delim("/workdir/mbb262/tf_binf_tu2021_v5.bed", header = F)
colnames(tf) <- c("seqid", "start", "end")

# add "chr" in front of all chromosome names (now its just numeric)
tf$seqid <- paste0("chr", tf$seqid)

# Turn into a genomic ranges object
gr_tf <- GenomicRanges::makeGRangesFromDataFrame(tf)

# Subset vcf file
tasGeno_tf <- tasGeno %>% 
  filterGenotypeTableSites(
    siteRangeFilterType = "none",
    gRangesObj = gr_tf)
# of SNPs: 845898

# calculate kinship
tasKin3 <- kinshipMatrix(tasObj = tasGeno_tf)


# 4 - K promoter + UTR ---------------------------------------------------------

# combine tf sites and 5 prime utr intervals
tf_strand <- tf # doesn't come with a strand, just say it's positive
tf_strand$strand <- rep("+", nrow(tf_strand))
tf_5utr <- rbind(tf_strand, gff_5utr)

# Turn into a genomic ranges object
gr_tf_5utr <- GenomicRanges::makeGRangesFromDataFrame(tf_5utr)

# Subset vcf file
tasGeno_tf_5utr <- tasGeno %>% 
  filterGenotypeTableSites(
    siteRangeFilterType = "none",
    gRangesObj = gr_tf_5utr)
# of SNPs: 871322

# calculate kinship
tasKin4 <- kinshipMatrix(tasObj = tasGeno_tf_5utr)


# 5 - K TEs within 2 kb of a gene ----------------------------------------------

# get gff file
system("wget -P /workdir/mbb262/ https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.TE.gff3.gz")

# Load in gff file, filter to 5' UTR regions of canonical transcripts
gff_te <- ape::read.gff("/workdir/mbb262/Zm-B73-REFERENCE-NAM-5.0.TE.gff3.gz", 
                        na.strings = c(".", "?"), GFF3 = TRUE) %>% 
  filter(seqid %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10")) %>% 
  select("seqid", "start", "end", "strand", "attributes")

# Take genes from gene gff, add 2kb to starts
gff_gene <- ape::read.gff("/workdir/mbb262/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz", 
                          na.strings = c(".", "?"), GFF3 = TRUE) %>% 
  filter(type == "gene", seqid %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10")) %>% 
  select("seqid", "start", "end", "strand", "attributes")
gff_gene$attributes <- gsub("ID=", "", gff_gene$attributes) # parse out gene names
gff_gene$attributes <- gsub(";.*", "", gff_gene$attributes) # parse out gene names

# Add a window
te_window <- 2000
gff_gene$start <- gff_gene$end - te_window

# Find overlapping TEs within a te_window size of a gene
gff_gene <- data.table::as.data.table(gff_gene)
gff_te <- data.table::as.data.table(gff_te)
setkey(gff_gene, seqid, start, end)
gene_te_overlaps <- data.table::foverlaps(gff_te, gff_gene, type="within",nomatch=0L) #0L drops TEs with no overlapping genes

# Select just the TE columns
gene_te_overlaps <- gene_te_overlaps %>% select("seqid", "i.start", "i.end", "i.strand", "attributes", "i.attributes")

# Fix column names
colnames(gene_te_overlaps) <- c("seqid", "start", "end", "strand", "v5_gene", "attributes")

# Turn into a genomic ranges object
gr_te <- GenomicRanges::makeGRangesFromDataFrame(gene_te_overlaps, keep.extra.columns = TRUE)

# Subset vcf file using Genomic Ranges object
tasGeno_te <- tasGeno %>% 
  filterGenotypeTableSites(
    siteRangeFilterType = "none",
    gRangesObj = gr_te)

# of sites left: 

# calculate kinship
tasKin5 <- kinshipMatrix(tasObj = tasGeno_te)


# Correlate matrices -----------------------------------------------------------

# Create matrices from tassel objects
tasKin1_matrix <- as.matrix(tasKin1)
tasKin2_matrix <- as.matrix(tasKin2)
tasKin3_matrix <- as.matrix(tasKin3)
tasKin4_matrix <- as.matrix(tasKin4)

# Needs to be re-written, create one table then correlate

# Correlate whole matrices to each other
cor(c(tasKin1_matrix), c(tasKin2_matrix)) # 0.9979598
cor(c(tasKin1_matrix), c(tasKin3_matrix)) # 0.9983401
cor(c(tasKin1_matrix), c(tasKin4_matrix)) # 0.9983531
cor(c(tasKin2_matrix), c(tasKin3_matrix)) # 0.9994387
cor(c(tasKin2_matrix), c(tasKin4_matrix)) # 0.9994798
cor(c(tasKin3_matrix), c(tasKin4_matrix)) # 0.9999985

# Gather lower off-diagonals 
# https://stackoverflow.com/questions/13049575/r-min-max-and-mean-of-off-diagonal-elements-in-a-matrix
kin1 <- tasKin1_matrix[row(tasKin1_matrix) == (col(tasKin1_matrix)+1)]
kin2 <- tasKin2_matrix[row(tasKin2_matrix) == (col(tasKin2_matrix)+1)]
kin3 <- tasKin3_matrix[row(tasKin3_matrix) == (col(tasKin3_matrix)+1)]
kin4 <- tasKin4_matrix[row(tasKin4_matrix) == (col(tasKin4_matrix)+1)]

# Correlate off-diagonals
cor(kin1, kin2) # 0.9989376
cor(kin1, kin3) # 0.9991454
cor(kin1, kin4) # 0.9991476
cor(kin2, kin3) # 0.9997494
cor(kin2, kin4) # 0.9997688
cor(kin3, kin4) # 0.9999993