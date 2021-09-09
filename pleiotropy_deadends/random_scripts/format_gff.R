# ---------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2020-07-09
#
# Description 
#   - Parse GFF file from Ensembl plants
# ---------------------------------------------------------------

# Load helpful packages
library(dplyr)

# Load in gff file 
# Genome gff dowloaded from here:
# ftp://ftp.ensemblgenomes.org/pub/plants/release-47/gff3/zea_mays/
maize_gff <- ape::read.gff("~/Downloads/Zea_mays.B73_RefGen_v4.47.gff3", na.strings = c(".", "?"), GFF3 = TRUE)

# Count total rows
total_rows <- nrow(maize_gff)
# Look at all levels in the gff file
levels(maize_gff$type)

# Parse gff - chromosomes
chroms <- maize_gff %>% filter(type == "chromosome")
chroms$transcript <- rep(NA, nrow(chroms)) # make dummy transcript column
chroms$gene <- rep(NA, nrow(chroms)) # make dummy gene column

# Parse gff - genes
genes <- maize_gff %>% filter(type %in% c("gene", "ncRNA_gene", "pseudogene"))
genes$transcript <- gsub("ID=gene:", "", genes$attributes) # extract genes, leave in column named transcript
genes$gene <- gsub(";.*", "", genes$transcript) # remove stuff after gene name
genes$transcript <- rep(NA, nrow(genes)) # replace transcripts column with NAs

# Create intergenic regions
gr_genes <- GenomicRanges::makeGRangesFromDataFrame(genes, keep.extra.columns = TRUE)
intergenic_region <- data.frame(GenomicRanges::gaps(gr_genes)) %>% select(-width, -strand)
intergenic_region$source <- rep("GRanges_gaps_between_gramene_genes", nrow(intergenic_region))
intergenic_region$type <- rep("intergenic_region", nrow(intergenic_region))
intergenic_region$intergenic_id <- paste0("intergenic_interval_", seq(1:nrow(intergenic_region)))
colnames(intergenic_region)[1] <- "seqid"

# Parse different RNAs and CDSs
rnas <- maize_gff %>% filter(type %in% c("CDS", "lnc_RNA", "miRNA", "mRNA", "ncRNA", "pre_miRNA", "pseudogenic_transcript", "RNase_MRP_RNA","rRNA", "snoRNA","snRNA" ,"SRP_RNA", "tRNA"))
rnas$transcript <- gsub("ID=transcript:", "", rnas$attributes) # extract transcripts
rnas$transcript <- gsub("ID=CDS:", "", rnas$transcript) # extract transcripts
rnas$transcript <- gsub(";.*", "", rnas$transcript) # extract transcript
rnas$gene <- gsub(".*Parent=gene:", "", rnas$attributes) # Get Parent gene
rnas$gene <- gsub(".*Parent=transcript:", "", rnas$gene) # Get Parent gene
rnas$gene <- gsub("_.*", "", rnas$gene) # Get Parent gene
rnas$gene <- gsub(";.*", "", rnas$gene) # get parent gene

# Parse exons, create introns, create a name column
exons <- maize_gff %>% filter(type == "exon")
exons$transcript <- gsub("Parent=transcript:", "", exons$attributes) # extract transcripts
exons$transcript <- gsub(";.*", "", exons$transcript) # extract transcript
exons$gene <- gsub("_T.*", "", exons$transcript) # Extract parent gene
exons$gene <- gsub("-T.*", "", exons$gene) # Extract parent gene
exons$exon_id <- gsub(".*Name=", "", exons$attributes) # extract transcripts
exons$exon_id <- gsub(";.*", "", exons$exon_id) # extract transcript

# Check if parsed correctly
exons[which.min(nchar(exons$gene)),10:12]
exons[which.max(nchar(exons$gene)),10:12]
exons[which.min(nchar(exons$transcript)),10:12]
exons[which.max(nchar(exons$transcript)),10:12]
exons[which.min(nchar(exons$exon_id)),10:12]
exons[which.max(nchar(exons$exon_id)),10:12]

# Parse UTRs
utrs <- maize_gff %>% filter(type %in% c("five_prime_UTR", "three_prime_UTR"))
utrs$transcript <- gsub("Parent=transcript:", "", utrs$attributes)
utrs$gene <- gsub("_T.*", "", utrs$transcript)

# rbind.fill all together and sort
together <- plyr::rbind.fill(chroms, genes, intergenic_region, rnas, exons, utrs) %>% arrange(seqid, start)

# Write to file
write.csv(together, "~/Downloads/Zea_mays.B73_RefGen_v4.47_parsed.csv", row.names = F, quote = F)

# Get only gene and intergene regions
gene_intergene <- together %>% filter(type %in% c("intergenic_region", "gene"))
write.csv(gene_intergene, "~/Downloads/Zea_mays.B73_RefGen_v4.47_gene_intergene_regions.csv",
          row.names = F, quote = F)

temp <- read.csv("~/Downloads/Zea_mays.B73_RefGen_v4.47_gene_intergene_regions.csv", header = T)


# -----------------------------
# Format open chromation data
# -----------------------------

# Load in Eli's 2015 open chromatin data
# Used Ensembl's Assembly converter to translate v3 coordinates to v4 coordinates
# https://plants.ensembl.org/Zea_mays/Tools/AssemblyConverter?db=core
shoot <- read.delim("~/Downloads/AP.bfthresh1.1.MNaseHS.Ranges_v4.bed", header = FALSE, sep = "\t")
root <- read.delim("~/Downloads/RP.bfthresh1.1.MNaseHS.Ranges_v4.bed", header = FALSE,sep = "\t")

# Add column names
colnames(shoot) <- c("chrom", "start", "end")
colnames(root) <- c("chrom", "start", "end")

# Make a dummy column and assign names --> named intervals needed for askDB
shoot$id <- paste0("shoot_open_chromatin_interval_", seq(1:nrow(shoot)))
root$id <- paste0("root_open_chromatin_interval_", seq(1:nrow(root)))

# Turn into a genomic ranges object
gr_shoot <- GenomicRanges::makeGRangesFromDataFrame(shoot, keep.extra.columns = TRUE)
gr_root <- GenomicRanges::makeGRangesFromDataFrame(root, keep.extra.columns = TRUE)

# Get all intervals not represented by the above intervals (think getting the intergenic regions to the genic)
# Also remove columns not needed
shoot_closed <- data.frame(GenomicRanges::gaps(gr_shoot)) %>% select(-width, -strand)
root_closed <- data.frame(GenomicRanges::gaps(gr_root)) %>% select(-width, -strand)

# Make a dummy column and fill with names
shoot_closed$id <- paste0("shoot_closed_chromatin_interval_", seq(1:nrow(shoot_closed)))
root_closed$id <- paste0("root_closed_chromatin_interval_", seq(1:nrow(root_closed)))

# Rename first column
colnames(shoot_closed)[1] <- "chrom"
colnames(root_closed)[1] <- "chrom"

# Rbind all together
shoot <- rbind(shoot, shoot_closed)
root <- rbind(root, root_closed)

# Save files
write.csv(shoot, "~/Downloads/AP.bfthresh1.1.MNaseHS.Ranges_v4_openclosedchromatin_intervals.csv", row.names = F, quote = F)
write.csv(root, "~/Downloads/RP.bfthresh1.1.MNaseHS.Ranges_v4_openclosedchromatin_intervals.csv", row.names = F, quote = F)



