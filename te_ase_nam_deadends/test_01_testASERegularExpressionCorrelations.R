# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-06-29
# Updated... 2023-06-30
#
# Description:
# Correlating expected vs observed expression counts
# ------------------------------------------------------------------------------

# # Load packages
# library(dplyr)
# library(ggplot2)
# 
# # Load in pan-gene table, melt
# pan <- read.csv("/workdir/mbb262/sim_reads/r_simulated/pan_gene_matrix_v3_cyverse.csv") %>% 
#   select("Pan_gene_id", "Rep_transcript", "B73", "Ky21", "CML247")
# pan_melt <- reshape2::melt(pan, id = c("Pan_gene_id", "Rep_transcript"))
# 
# # Read in files
# b73_sim_counts_matrix_expected <- read.csv("/workdir/mbb262/sim_reads/r_simulated/b73_sim_counts.csv")
# ky21_sim_counts_matrix_expected <- read.csv("/workdir/mbb262/sim_reads/r_simulated/ky21_sim_counts.csv")
# cml247_sim_counts_matrix_expected <- read.csv("/workdir/mbb262/sim_reads/r_simulated/cml247_sim_counts.csv")
# hybridB73Ky21_sim_counts_matrix_expected <- read.csv("/workdir/mbb262/sim_reads/r_simulated/hybridB73Ky21_sim_counts.csv")
# hybridB73Cml247_sim_counts_matrix_expected <- read.csv("/workdir/mbb262/sim_reads/r_simulated/hybridB73Cml247_sim_counts.csv")
# hybridCml247Ky21_sim_counts_matrix_expected <- read.csv("/workdir/mbb262/sim_reads/r_simulated/hybridCml247Ky21_sim_counts.csv")
# 
# # Load in observed counts
# hybridB73Ky21_observed <- read.delim("/workdir/mbb262/sim_reads/output/read_counts/hybridB73Ky21_minimap2.txt")
# hybridB73CML247_observed <- read.delim("/workdir/mbb262/sim_reads/output/read_counts/hybridB73Cml247_minimap2.txt")
# hybridCML247Ky21_observed <- read.delim("/workdir/mbb262/sim_reads/output/read_counts/hybridCml247Ky21_minimap2.txt")
# inbred_observed <- read.delim("/workdir/mbb262/sim_reads/output/read_counts/inbred_counts_minimap2.txt")
# 
# # Format inbred counts
# trans_inbred_observed <- t(inbred_observed[,-1]) %>% data.frame()
# colnames(trans_inbred_observed) <- inbred_observed[,1]
# trans_inbred_observed$transcript <- rownames(trans_inbred_observed)
# rownames(trans_inbred_observed) <- NULL
# trans_inbred_observed$gene <- gsub("_.*", "", trans_inbred_observed$transcript)
# trans_inbred_observed <- trans_inbred_observed %>% 
#   select("gene", "transcript", "b73_to_b73_minimap2", "ky21_to_ky21_minimap2", "cml247_to_cml247_minimap2")
# 
# 
# # Merge expected counts with observed ------------------------------------------
# 
# # Inbred B73 --------------------------------
# inbred_b73 <- merge(b73_sim_counts_matrix_expected, trans_inbred_observed %>% select("gene", "b73_to_b73_minimap2"), 
#                     by = "gene")
# cor(inbred_b73$sample_01, inbred_b73$b73_to_b73_minimap2, method = "spearman")^2
# a <- ggplot(inbred_b73, aes(x = sample_01, y = b73_to_b73_minimap2)) +
#   geom_point() +
#   xlab("Expected (aka Simulated)") +
#   ylab("Observed minimap Kotlin") +
#   ggtitle("Inbred B73") +
#   geom_abline(slope = 1)+
#   geom_text(x= 15, y = 65, label = paste0("R2 = ", round(cor(inbred_b73$sample_01, inbred_b73$b73_to_b73_minimap2, method = "spearman")^2,4)))
# 
# # Inbred Ky21 ------------------------------
# pangene_ky21 <- merge(pan, trans_inbred_observed %>% select("transcript", "ky21_to_ky21_minimap2"), 
#                       by.x = "B73", by.y = "transcript")
# dim(pangene_ky21)
# 
# pangene_ky21$cleangene <- gsub("_.*", "", pangene_ky21$Ky21)
# inbred_ky21 <- merge(ky21_sim_counts_matrix_expected, pangene_ky21, 
#                      by.x = "gene", by.y = "cleangene")
# 
# cor(inbred_ky21$sample_01, inbred_ky21$ky21_to_ky21_minimap2, method = "spearman")^2
# b <- ggplot(inbred_ky21, aes(x = sample_01, y = ky21_to_ky21_minimap2)) +
#   geom_point() +
#   xlab("Expected (aka Simulated)") +
#   ylab("Observed minimap2 Kotlin") +
#   ggtitle("Inbred Ky21") +
#   geom_abline(slope = 1)+
#   geom_text(x= 15, y =150, label = paste0("R2 = ", round(cor(inbred_ky21$sample_01, inbred_ky21$ky21_to_ky21_minimap2, method = "spearman")^2,4)))
# 
# 
# # Inbred CML247 ------------------------------
# pangene_cml247 <- merge(pan, trans_inbred_observed %>% select("transcript", "cml247_to_cml247_minimap2"), 
#                         by.x = "B73", by.y = "transcript")
# dim(pangene_cml247)
# 
# pangene_cml247$cleangene <- gsub("_.*", "", pangene_cml247$CML247)
# inbred_cml247 <- merge(cml247_sim_counts_matrix_expected, pangene_cml247, 
#                        by.x = "gene", by.y = "CML247")
# 
# cor(inbred_cml247$sample_01, inbred_cml247$cml247_to_cml247_minimap2, method = "spearman")^2
# c <- ggplot(inbred_cml247, aes(x = sample_01, y = cml247_to_cml247_minimap2)) +
#   geom_point() +
#   xlab("Expected (aka Simulated)") +
#   ylab("Observed minimap2 Kotlin") +
#   ggtitle("Inbred CML247") +
#   geom_abline(slope = 1)+
#   geom_text(x= 15, y =110, label = paste0("R2 = ", round(cor(inbred_cml247$sample_01, inbred_cml247$cml247_to_cml247_minimap2, method = "spearman")^2,4)))
# 
# library(patchwork)
# (a+b)/c
# # ggsave("/home/mbb262/git_projects/te_ase_nam/figs/inbred_counter.png", plot = (a+b)/(c+plot_spacer()), width = 9, height = 5)
# 
# # Hybrid B73xKy21 to B73 ----------------------
# hbk2b_ky21 <- merge(pan_melt, hybridB73Ky21_observed, 
#                     by.x = "value", by.y = "transcript_id")
# 
# 
# 
# 
# # -------------------------------------------------
# # Correlate salmon with expected
# # ------------------------------------------------
# 
# # Load packages
# library(dplyr)
# library(ggplot2)
# 
# # Load in pan-gene table, melt
# pan <- read.csv("/workdir/mbb262/sim_reads/r_simulated/pan_gene_matrix_v3_cyverse.csv") %>% 
#   select("Pan_gene_id", "Rep_transcript", "B73", "Ky21", "CML247")
# pan_melt <- reshape2::melt(pan, id = c("Pan_gene_id", "Rep_transcript"))
# 
# # Read in files
# b73_sim_counts_matrix_expected <- read.csv("/workdir/mbb262/sim_reads/r_simulated/b73_sim_counts.csv")
# ky21_sim_counts_matrix_expected <- read.csv("/workdir/mbb262/sim_reads/r_simulated/ky21_sim_counts.csv")
# cml247_sim_counts_matrix_expected <- read.csv("/workdir/mbb262/sim_reads/r_simulated/cml247_sim_counts.csv")
# # hybridB73Ky21_sim_counts_matrix_expected <- read.csv("/workdir/mbb262/sim_reads/r_simulated/hybridB73Ky21_sim_counts.csv")
# # hybridB73Cml247_sim_counts_matrix_expected <- read.csv("/workdir/mbb262/sim_reads/r_simulated/hybridB73Cml247_sim_counts.csv")
# # hybridCml247Ky21_sim_counts_matrix_expected <- read.csv("/workdir/mbb262/sim_reads/r_simulated/hybridCml247Ky21_sim_counts.csv")
# 
# # read in salmon counts
# b73_observed <- read.delim("/workdir/mbb262/sim_reads/salmon/b73/b73_transcripts_quant/quant.sf")
# b73_observed$Name <- gsub("_.*", "", b73_observed$Name)
# 
# cml247_observed <- read.delim("/workdir/mbb262/sim_reads/salmon/cml247/cml247_transcripts_quant/quant.sf")
# cml247_observed$Name <- gsub("_.*", "", cml247_observed$Name)
# 
# ky21_observed <- read.delim("/workdir/mbb262/sim_reads/salmon/ky21/ky21_transcripts_quant/quant.sf")
# ky21_observed$Name <- gsub("_.*", "", ky21_observed$Name)
# 
# 
# # Merge and correlate b73 ---------------------------
# temp <- merge(b73_sim_counts_matrix_expected, b73_observed, by.x = "gene", by.y = "Name")
# cor(temp$sample_01, temp$NumReads, method = "spearman")^2
# a <- ggplot(temp, aes(x = sample_01, y = NumReads)) +
#   geom_point() +
#   xlab("Expected (aka Simulated)") +
#   ylab("Observed salmon") +
#   ggtitle("Inbred B73") +
#   geom_abline(slope = 1)+
#   geom_text(x= 15, y = 170, label = paste0("R2 = ", round(cor(temp$sample_01, temp$NumReads, method = "spearman")^2,4)))
# 
# # ggsave("/home/mbb262/git_projects/te_ase_nam/figs/b73_salmon.png", plot = a, width = 7, height = 5)
# 
# 
# # Inbred Ky21 ------------------------------
# pangene_ky21 <- merge(pan, ky21_observed, 
#                       by.x = "Pan_gene_id", by.y = "Name")
# dim(pangene_ky21)
# 
# pangene_ky21$cleangene <- gsub("_.*", "", pangene_ky21$Ky21)
# inbred_ky21 <- merge(ky21_sim_counts_matrix_expected, pangene_ky21, 
#                      by.x = "gene", by.y = "cleangene")
# 
# cor(inbred_ky21$sample_01, inbred_ky21$NumReads, method = "spearman")^2
# b <- ggplot(inbred_ky21, aes(x = sample_01, y = NumReads)) +
#   geom_point() +
#   xlab("Expected (aka Simulated)") +
#   ylab("Observed salmon") +
#   ggtitle("Inbred Ky21") +
#   geom_abline(slope = 1)+
#   geom_text(x= 15, y =170, label = paste0("R2 = ", round(cor(inbred_ky21$sample_01, inbred_ky21$NumReads, method = "spearman")^2,4)))
# 
# 
# # Inbred CML247 ------------------------------
# pangene_cml247 <- merge(pan, cml247_observed, 
#                         by.x = "Pan_gene_id", by.y = "Name")
# dim(pangene_cml247)
# 
# pangene_cml247$cleangene <- gsub("_.*", "", pangene_cml247$CML247)
# inbred_cml247 <- merge(cml247_sim_counts_matrix_expected, pangene_cml247, 
#                        by.x = "gene", by.y = "CML247")
# 
# cor(inbred_cml247$sample_01, inbred_cml247$NumReads, method = "spearman")^2
# c <- ggplot(inbred_cml247, aes(x = sample_01, y = NumReads)) +
#   geom_point() +
#   xlab("Expected (aka Simulated)") +
#   ylab("Observed salmon") +
#   ggtitle("Inbred CML247") +
#   geom_abline(slope = 1)+
#   geom_text(x= 15, y =110, label = paste0("R2 = ", round(cor(inbred_cml247$sample_01, inbred_cml247$NumReads, method = "spearman")^2,4)))
# 
# 
# # Save plot
# library(patchwork)
# (a+b)/c
# # ggsave("/home/mbb262/git_projects/te_ase_nam/figs/inbred_counter.png", plot = (a+b)/(c+plot_spacer()), width = 9, height = 5)


# -------------------------------------------------------------------------------
# Test SEESAW
# https://mikelove.github.io/Bioc2022AllelicExpression/articles/Bioc2022AllelicExpression.html
# ------------------------------------------------------------------------------

# Load packages
library(GenomicRanges)
library(fishpond)
library(dplyr)

# Filter quantified transcripts to only those found between both genomes
temp <- read.delim("/workdir/mbb262/sim_reads/salmon/hybridb73ky21//hybrid_transcripts_quant/quant.sf")
temp$gene <- gsub("_ky21", "", temp$Name)
temp$gene <- gsub("_b73", "", temp$gene)
counts_gene <- table(temp$gene) %>% data.frame() %>% filter(Freq == 2) %>% select(Var1)
counts_gene$ky21 <- paste0(counts_gene$Var1, "_ky21")
counts_gene$b73 <- paste0(counts_gene$Var1, "_b73")
all_genes_ids_to_keep <- c(as.vector(counts_gene$b73), as.vector(counts_gene$ky21))
temp <- temp %>% filter(Name %in% all_genes_ids_to_keep) %>% select(-gene)
write.table(temp, 
            "/workdir/mbb262/sim_reads/salmon/hybridb73ky21/hybrid_transcripts_quant/shared_quant.sf",
            quote = F, row.names = F, sep = '\t')


# Get coordinates for transcripts
# wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz
b73_gff <- ape::read.gff("/workdir/mbb262/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3") %>% 
  dplyr::filter(type == "mRNA" & seqid %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10")) %>% 
  dplyr::select("seqid", "start", "end", "strand", "attributes")
b73_gff$tx_id <- gsub("ID=", "", b73_gff$attributes)
b73_gff$tx_id <- gsub(";.*", "", b73_gff$tx_id)
b73_gff$gene_id <- gsub(".*Parent=", "", b73_gff$attributes)
b73_gff$gene_id <- gsub(";.*", "", b73_gff$gene_id)
b73_gff <- b73_gff %>% dplyr::select(-"attributes")

# Turn into a genomic ranges object
combined_gff_gr <- GenomicRanges::makeGRangesFromDataFrame(b73_gff, keep.extra.columns = TRUE)



# # Debugging source code ------------------------------------------------------
# temp <- read.delim("/workdir/mbb262/sim_reads/salmon/hybridb73ky21//hybrid_transcripts_quant/shared_quant.sf")
# table(grepl("_b73", temp$Name)) # True = B73, False = other genome
# 
# txOut <- is.null(t2g)
# se <- tximeta::tximeta("/workdir/mbb262/sim_reads/salmon/hybridb73ky21/hybrid_transcripts_quant/shared_quant.sf",
#                        skipMeta=TRUE, tx2gene=t2g)
# 
# # remove any characters after "|"
# rownames(se) <- sub("\\|.*", "", rownames(se))
# 
# ntxp <- nrow(se)/2
# n <- ncol(se)
# 
# # gather transcript names for a1 and a2 alleles
# a1match <- paste0("_","ky21","$")
# a2match <- paste0("_","b73","$")
# 
# txp_nms_a1 <- grep(a1match, rownames(se), value=TRUE)
# stopifnot(length(txp_nms_a1) == ntxp)

# Finish debugging source code -------------------------------------------------

# Link transcripts with genes (tx2gene)
t2g <- fishpond::makeTx2Tss(combined_gff_gr) # GRanges object
mcols(t2g)[,c("tx_id","group_id")]

# `files` points to `quant.sf` files in original Salmon directories
files <- fishpond::importAllelicCounts("/workdir/mbb262/sim_reads/salmon/hybridb73ky21/hybrid_transcripts_quant/shared_quant.sf",
                                       a1 = 'ky21', a2 = 'b73', tx2gene = t2g, dropInfReps = TRUE)

SummarizedExperiment::colData(files)
files <- labelKeep(files)
files <- swish(files, x="allele", pair="sample", cor="pearson")

# Differential AI 
set.seed(1)
y <- makeSimSwishData(allelic = T, dynamic = F, n=12)
colData(y)
table(y$condition, y$allele)
y <- labelKeep(y)
y <- swish(y, x="allele", pair="sample",
           cov=NULL, interaction=TRUE)
mcols(y)[1:2,c("stat","qvalue")]


# Simulate data for dynamic allelic imbalance section
set.seed(1)
y <- makeSimSwishData(dynamic=TRUE)
SummarizedExperiment::colData(y)
y <- labelKeep(y)
y <- swish(y, x="allele", pair="sample", cov="time", cor="pearson")
mcols(y)[1:2,c("stat","qvalue")]








