# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-11-22
# Updated... 2023-11-22
#
# Description:
# Test pipeline to subset fastas to the same set of ids before mapping with salmon
# ------------------------------------------------------------------------------

# R: Read all files into R, make into a list
cml322 <- read.delim("/workdir/mbb262/test_intersect_transcripts/cml322.txt", header = F)
phb47 <- read.delim("/workdir/mbb262/test_intersect_transcripts/phb47.txt", header = F)

# R: Intersect transcript names
lala <- data.frame(intersect(cml322$V1, phb47$V1))

# R: Check if lengths changed
dim(cml322)
dim(phb47)
dim(lala)

# R: Run another check to make sure they're present in more than 1 genome
temp <- rbind(cml322, phb47)
countTrans <- table(temp$V1)
countLala <- table(lala) # should all = 1

# R: Export common set of ids
write.table(lala, "/workdir/mbb262/test_intersect_transcripts/cmlphbshared.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE)
    

# Try loading this data into R -------------------------------------------------

# Load packages
library(dplyr)


## Key error where this all begins ---------------------------------------------

# Find files
samples <- list.files("/workdir/mbb262/test_intersect_transcripts/output", full.names = TRUE)
files <- file.path(samples, "quant.sf")

# Subset to debug --> c("MS21R193", "MS21R254", "MS21R317", "MS21R384")
samples <- samples[c(60,120, 183, 248)]
files <- files[c(60,120, 183, 248)]

names(files) <- gsub(".*/", "", samples)
files_df <- data.frame(files)
files_df$sample_name <- rownames(files_df)

# Make metadata for loading into R
meta_hybrids <- read.csv("/home/mbb262/git_projects/te_ase_nam/data/sample_metadata.csv") %>% 
  dplyr::filter(type == "hybrid") %>% 
  dplyr::select(sample_name, cultivar, group, age, tissue)
temp <- table(meta_hybrids$cultivar)
coldata <- merge(meta_hybrids, files_df, by = "sample_name") 
colnames(coldata) <- c("names", "cultivar", "group", "age", "tissue", "files")
dim(coldata)
head(coldata)

# Load gff, make transcript file, turn into genomic ranges object, then makeTx2Tss object
b73_gff <- ape::read.gff("/workdir/mbb262/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3") %>% 
  dplyr::filter(type == "mRNA" & seqid %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10")) %>% 
  dplyr::select("seqid", "start", "end", "strand", "attributes")
b73_gff$tx_id <- gsub("ID=", "", b73_gff$attributes)
b73_gff$tx_id <- gsub(";.*", "", b73_gff$tx_id)
b73_gff$geneID <- gsub(".*Parent=", "", b73_gff$attributes)
b73_gff$geneID <- gsub(";.*", "", b73_gff$geneID)
b73_gff <- b73_gff %>% dplyr::select(-"attributes")
colnames(b73_gff)[6] <- "gene_id"
combined_gff_gr <- GenomicRanges::makeGRangesFromDataFrame(b73_gff, keep.extra.columns = TRUE)
t2g <- fishpond::makeTx2Tss(combined_gff_gr) # GRanges object
GenomicRanges::mcols(t2g)[,c("tx_id","group_id")]

# Try loading in with seesaw
y <- fishpond::importAllelicCounts(coldata,
                                   a1="PHB47", 
                                   a2="CML322-REFERENCE-NAM-1.0_agat_sub.fa",
                                   format="assays", 
                                   tx2gene=t2g,
                                   dropInfReps = FALSE)



# For all files ----------------------------------------------------------------

# R: Read all files into R, make into a list
all_files <- list.files(path = "/workdir/mbb262/nam_hybrid_rnaseq/salmon/nam_genomes/liftoff_inbred/agat_sequences/agat_transcript_names",
                        pattern = ".txt", full.names = TRUE)
all_files_list <- lapply(all_files, read.delim, header = F)
all_files_list <- lapply(all_files_list, function (x) {as.vector(x$V1)})

# R: Intersect transcript names
lala <- data.frame(Reduce(intersect, all_files_list))
head(lala)

# R: Check if lengths changed
lapply(all_files_list, length)
dim(lala)

# R: Export common set of ids
write.table(lala, "/workdir/mbb262/nam_hybrid_rnaseq/salmon/nam_genomes/liftoff_inbred/agat_sequences/agat_transcript_names/shared_transcripts_across_all_NAM.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE)




