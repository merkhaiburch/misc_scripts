# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-02-07
# Updated... 2023-02-13
#
# Description:
# Simulate reads from B73, Cml247, and an artificial B73xCml247 hybrid and
#                from B73, Ky21, and an artificial B73xKy21 hybrid 
# ------------------------------------------------------------------------------

# Load packages
library(dplyr)

# Grab only the transcript sequences from B73, Cml247 and Ky21
b73_gff <- ape::read.gff("/workdir/mbb262/sim_reads/genomes_to_simulate_reads_from/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz", 
                          na.strings = c(".", "?"), GFF3 = TRUE) %>% 
  filter(type == "mRNA", seqid %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10")) %>% 
  select("seqid", "start", "end", "strand", "attributes")

cml247_gff <- ape::read.gff("/workdir/mbb262/sim_reads/genomes_to_simulate_reads_from/Zm-CML247-REFERENCE-NAM-1.0_Zm00023ab.1.gff3.gz", 
                     na.strings = c(".", "?"), GFF3 = TRUE) %>% 
  filter(type == "mRNA", seqid %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10")) %>% 
  select("seqid", "start", "end", "strand", "attributes")

ky21_gff <- ape::read.gff("/workdir/mbb262/sim_reads/genomes_to_simulate_reads_from/Zm-Ky21-REFERENCE-NAM-1.0_Zm00031ab.1.gff3.gz", 
                          na.strings = c(".", "?"), GFF3 = TRUE) %>% 
  filter(type == "CDS", seqid %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10")) %>% 
  select("seqid", "start", "end", "strand", "attributes")

# Format attributes
b73_gff$gene <- gsub("ID=", "", b73_gff$attributes)
b73_gff$gene <- gsub(";.*", "", b73_gff$gene)
cml247_gff$gene <- gsub("ID=", "", cml247_gff$attributes)
cml247_gff$gene <- gsub(";.*", "", cml247_gff$gene)
ky21_gff$gene <- gsub(".*Parent=", "", ky21_gff$attributes)
ky21_gff$gene <- gsub(";.*", "", ky21_gff$gene)

# Remove old column 
b73_gff <- b73_gff %>% select(-"attributes")
cml247_gff <- cml247_gff %>% select(-"attributes")
ky21_gff <- ky21_gff %>% select(-"attributes")

# Group by gene, select the smallest start and largest end
cml247_transcript_start <- cml247_gff %>% group_by(gene) %>% slice(which.min(start)) %>% select(-"end")
cml247_transcript_end <- cml247_gff %>% group_by(gene) %>% slice(which.max(end)) %>% select(-"start")
ky21_transcript_start <- ky21_gff %>% group_by(gene) %>% slice(which.min(start)) %>% select(-"end")
ky21_transcript_end <- ky21_gff %>% group_by(gene) %>% slice(which.max(end)) %>% select(-"start")

# Merge together
cml247_transcript <- merge(cml247_transcript_start, cml247_transcript_end, by = c("gene", "strand", "seqid"))
ky21_transcript <- merge(ky21_transcript_start, ky21_transcript_end, by = c("gene", "strand", "seqid"))

# Select only canonical transcripts
ky21_canonical <- read.csv("/workdir/mbb262/sim_reads/genomes_to_simulate_reads_from/Zm-Ky21-REFERENCE-NAM-1.0_Zm00031ab.1.canonical_transcripts", header = F)
ky21_transcript <- ky21_transcript %>% filter(gene %in% ky21_canonical$V1)
max(table(unique(gsub("_T.*", "", ky21_transcript$gene)))) # subsetted to canonical, should = 1
dim(ky21_transcript)

b73_reel <- read.csv("/workdir/mbb262/sim_reads/genomes_to_simulate_reads_from/top_reel_transcripts.txt", header = F)
b73_transcript <- b73_gff %>% filter(gene %in% b73_reel$V1)
dim(b73_transcript)
max(table(unique(gsub("_T.*", "", b73_transcript$gene)))) # subsetted to reel, should = 1

cml247_canonical <- read.csv("/workdir/mbb262/sim_reads/genomes_to_simulate_reads_from/Zm-CML247-REFERENCE-NAM-1.0_Zm00023ab.1.canonical_transcripts", header = F)
cml247_transcript <- cml247_transcript %>% filter(gene %in% cml247_canonical$V1)
max(table(unique(gsub("_T.*", "", cml247_transcript$gene)))) # subsetted to canonical, should = 1
dim(cml247_transcript)


# Rearrange
b73_transcript <- b73_transcript %>% select("seqid", "start", "end", "strand", "gene")
cml247_transcript <- cml247_transcript %>% select("seqid", "start", "end", "strand", "gene")
ky21_transcript <- ky21_transcript %>% select("seqid", "start", "end", "strand", "gene")

# Export bed file
library(GenomicRanges)
b73_gr <- GenomicRanges::makeGRangesFromDataFrame(df = b73_transcript, keep.extra.columns = TRUE)
b73_grdf <- data.frame(seqnames=seqnames(b73_gr),
                        starts=start(b73_gr)-1,
                        ends=end(b73_gr),
                        name=b73_gr$gene)

write.table(b73_grdf, file="/workdir/mbb262/sim_reads/genomes_to_simulate_reads_from/b73_transcript.bed", 
            quote=F, sep="\t", row.names=F, col.names=F)

cml247_gr <- GenomicRanges::makeGRangesFromDataFrame(df = cml247_transcript, keep.extra.columns = TRUE)
cml247_grdf <- data.frame(seqnames=seqnames(cml247_gr),
                 starts=start(cml247_gr)-1,
                 ends=end(cml247_gr),
                 name=cml247_gr$gene)

write.table(cml247_grdf, file="/workdir/mbb262/sim_reads/genomes_to_simulate_reads_from/cml247_transcript.bed", 
            quote=F, sep="\t", row.names=F, col.names=F)

ky21_gr <- GenomicRanges::makeGRangesFromDataFrame(df = ky21_transcript, keep.extra.columns = TRUE)
ky21_grdf <- data.frame(seqnames=seqnames(ky21_gr),
                        starts=start(ky21_gr)-1,
                        ends=end(ky21_gr),
                        name=ky21_gr$gene)
write.table(ky21_grdf, file="/workdir/mbb262/sim_reads/genomes_to_simulate_reads_from/ky21_transcript.bed", 
            quote=F, sep="\t", row.names=F, col.names=F)


## Simulate reads --------------------------------------------------------------

# Load packages
library(polyester)
library(Biostrings)

# Set in and out paths
fastapath_s1 <- "/workdir/mbb262/sim_reads/genomes_to_simulate_reads_from/b73_transcripts_named.fasta"
fastapath_s2 <- "/workdir/mbb262/sim_reads/genomes_to_simulate_reads_from/cml247_transcripts_named.fasta"
fastapath_s3 <- "/workdir/mbb262/sim_reads/genomes_to_simulate_reads_from/ky21_transcripts_named.fasta"
fastapath_s4 <- "/workdir/mbb262/sim_reads/genomes_to_simulate_reads_from/b73_ky21_hybrid.fasta"
fastapath_s5 <- "/workdir/mbb262/sim_reads/genomes_to_simulate_reads_from/b73_cml247_hybrid.fasta"
fastapath_s6 <- "/workdir/mbb262/sim_reads/genomes_to_simulate_reads_from/cml247_ky21_hybrid.fasta"

outdir <- "/workdir/mbb262/sim_reads/r_simulated/"

# Count transcripts for all files in directory
s1_numtx <- polyester::count_transcripts(fastapath_s1)
s2_numtx <- polyester::count_transcripts(fastapath_s2)
s3_numtx <- polyester::count_transcripts(fastapath_s3)
s4_numtx <- polyester::count_transcripts(fastapath_s4)
s5_numtx <- polyester::count_transcripts(fastapath_s5)
s6_numtx <- polyester::count_transcripts(fastapath_s6)

# Format names
s1_tNames <-  gsub("'", "", names(Biostrings::readDNAStringSet(fastapath_s1)))
s2_tNames <-  gsub("'", "", names(Biostrings::readDNAStringSet(fastapath_s2)))
s3_tNames <-  gsub("'", "", names(Biostrings::readDNAStringSet(fastapath_s3)))
s4_tNames <-  gsub("'", "", names(Biostrings::readDNAStringSet(fastapath_s4)))
s5_tNames <-  gsub("'", "", names(Biostrings::readDNAStringSet(fastapath_s5)))
s6_tNames <-  gsub("'", "", names(Biostrings::readDNAStringSet(fastapath_s6)))

# Simulate experiments
s1_res <- polyester::simulate_experiment(fastapath_s1, 
                           reads_per_transcript=20,
                           fold_changes=rep(1, s1_numtx), 
                           outdir= paste0(outdir, "/b73"),
                           transcriptid=s1_tNames, 
                           num_reps = 1,
                           readlen = 150,
                           paired = FALSE)

s2_res <- polyester::simulate_experiment(fastapath_s2, 
                                         reads_per_transcript=20,
                                         fold_changes=rep(1, s2_numtx), 
                                         outdir= paste0(outdir, "/cml247"),
                                         transcriptid=s2_tNames, 
                                         num_reps = 1,
                                         readlen = 150,
                                         paired = FALSE)

s3_res <- polyester::simulate_experiment(fastapath_s3, 
                                         reads_per_transcript=20,
                                         fold_changes=rep(1, s3_numtx), 
                                         outdir= paste0(outdir, "/ky21"),
                                         transcriptid=s3_tNames, 
                                         num_reps = 1,
                                         readlen = 150,
                                         paired = FALSE)

s4_res <- polyester::simulate_experiment(fastapath_s4, 
                                         reads_per_transcript=20,
                                         fold_changes=rep(1, s4_numtx), 
                                         outdir= paste0(outdir, "/hybridB73Cml247"),
                                         transcriptid=s4_tNames, 
                                         num_reps = 1,
                                         readlen = 150,
                                         paired = FALSE)

s5_res <- polyester::simulate_experiment(fastapath_s5, 
                                         reads_per_transcript=20,
                                         fold_changes=rep(1, s5_numtx), 
                                         outdir= paste0(outdir, "/hybridB73Ky21"),
                                         transcriptid=s5_tNames, 
                                         num_reps = 1,
                                         readlen = 150,
                                         paired = FALSE)

s6_res <- polyester::simulate_experiment(fastapath_s6, 
                                         reads_per_transcript=20,
                                         fold_changes=rep(1, s6_numtx), 
                                         outdir= paste0(outdir, "/hybridCml247Ky21"),
                                         transcriptid=s6_tNames, 
                                         num_reps = 1,
                                         readlen = 150,
                                         paired = FALSE)

# rename files in each directory because polyester wont do it for me


## Turn .rds files into dataframes for easy correlations -----------------------

library(dplyr)

load("/workdir/mbb262/sim_reads/r_simulated/b73/sample_01_b73.rda")
b73_sim_counts_matrix <- counts_matrix %>% data.frame()
b73_sim_counts_matrix$gene <- gsub("_.*", "", rownames(b73_sim_counts_matrix))
head(b73_sim_counts_matrix)
dim(b73_sim_counts_matrix)

load("/workdir/mbb262/sim_reads/r_simulated/ky21/sample_01_ky21.rda")
ky21_sim_counts_matrix <- counts_matrix %>% data.frame()
ky21_sim_counts_matrix$gene <- gsub("_.*", "", rownames(ky21_sim_counts_matrix))
head(ky21_sim_counts_matrix)
dim(ky21_sim_counts_matrix)

load("/workdir/mbb262/sim_reads/r_simulated/cml247/sample_01_cml247.rda")
cml247_sim_counts_matrix <- counts_matrix %>% data.frame()
cml247_sim_counts_matrix$gene <- gsub("_cml247", "", rownames(cml247_sim_counts_matrix))
head(cml247_sim_counts_matrix)
dim(cml247_sim_counts_matrix)

load("/workdir/mbb262/sim_reads/r_simulated/hybridB73Ky21/sample_01_hybridB73Ky21.rda")
hybridB73Ky21_sim_counts_matrix <- counts_matrix %>% data.frame()
hybridB73Ky21_sim_counts_matrix$gene <- gsub("_.*", "", rownames(hybridB73Ky21_sim_counts_matrix))
head(hybridB73Ky21_sim_counts_matrix)
dim(hybridB73Ky21_sim_counts_matrix)

load("/workdir/mbb262/sim_reads/r_simulated/hybridB73Cml247/sample_01_hybridB73Cml247.rda")
hybridB73Cml247_sim_counts_matrix <- counts_matrix %>% data.frame()
hybridB73Cml247_sim_counts_matrix$gene <- gsub("_cml247", "", rownames(hybridB73Cml247_sim_counts_matrix))
hybridB73Cml247_sim_counts_matrix$gene <- gsub("_b73", "", hybridB73Cml247_sim_counts_matrix$gene)
head(hybridB73Cml247_sim_counts_matrix)
dim(hybridB73Cml247_sim_counts_matrix)

load("/workdir/mbb262/sim_reads/r_simulated/hybridCml247Ky21/sample_01_hybridCml247Ky21.rda")
hybridCml247Ky21_sim_counts_matrix <- counts_matrix %>% data.frame()
hybridCml247Ky21_sim_counts_matrix$gene <- gsub("_cml247", "", rownames(hybridCml247Ky21_sim_counts_matrix))
hybridCml247Ky21_sim_counts_matrix$gene <- gsub("_ky21", "", hybridCml247Ky21_sim_counts_matrix$gene)
head(hybridCml247Ky21_sim_counts_matrix)
dim(hybridCml247Ky21_sim_counts_matrix)

# Export to file
write.csv(b73_sim_counts_matrix, "/workdir/mbb262/sim_reads/r_simulated/b73_sim_counts.csv",
          row.names = FALSE, quote = FALSE)

write.csv(ky21_sim_counts_matrix, "/workdir/mbb262/sim_reads/r_simulated/ky21_sim_counts.csv",
          row.names = FALSE, quote = FALSE)

write.csv(cml247_sim_counts_matrix, "/workdir/mbb262/sim_reads/r_simulated/cml247_sim_counts.csv",
          row.names = FALSE, quote = FALSE)

write.csv(hybridB73Ky21_sim_counts_matrix, "/workdir/mbb262/sim_reads/r_simulated/hybridB73Ky21_sim_counts.csv",
          row.names = FALSE, quote = FALSE)

write.csv(hybridB73Cml247_sim_counts_matrix, "/workdir/mbb262/sim_reads/r_simulated/hybridB73Cml247_sim_counts.csv",
          row.names = FALSE, quote = FALSE)

write.csv(hybridCml247Ky21_sim_counts_matrix, "/workdir/mbb262/sim_reads/r_simulated/hybridCml247Ky21_sim_counts.csv",
          row.names = FALSE, quote = FALSE)

# Make a file with all transcript names for counting
trans_ids <- c(gsub("_b73", "", rownames(b73_sim_counts_matrix)),
                   gsub("_ky21", "", rownames(ky21_sim_counts_matrix)),
                   gsub("_cml247", "", rownames(cml247_sim_counts_matrix))) %>% 
  data.frame()
colnames(trans_ids) <- NULL
write.csv(trans_ids, "/workdir/mbb262/sim_reads/r_simulated/simulated_transcript_ids.txt",
          row.names = FALSE, quote = FALSE)






