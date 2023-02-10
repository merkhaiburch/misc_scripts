# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2021-12-21
#
# Description 
#   - 06_count_normalization.R
#   - Normalize feature count data from HTSeq
# Detailed Purpose:
#    The main purpose of this Rscript is to normalize features that
#    have been generated with HTSeq. Additionally, this script will
#    be used to filter the matrix based on bad taxa IDs (see `07`)
#    and gene IDs with low counts (<= 10 counts) across all taxa
#    samples.
#
# NOTE:
#    Scripts `06` and `07` work in conjunction.
#--------------------------------------------------------------------

# === Preamble ======================================================

## Load packages ----
library(data.table)
library(DESeq2)
library(dplyr)
library(magrittr)
library(pbapply)
library(readr)


## Parameters ----
count_path <- "/workdir/mbb262/nam_hybrid_rnaseq/output/htseq_counts/"


## Load data ----
# NOTE: data is loaded via DESeq2 - much simpler for downstream...



# === Process data ==================================================

## Get Taxa IDs ----
taxa_ids <- list.files(
  path = count_path,
  pattern = "\\.counts$",
  full.names = TRUE
) %>%
  stringr::str_replace(pattern = "^.*/.*/", replacement = "") %>%
  stringr::str_replace(pattern = "_CKD.*$", replacement = "")



# === DESeq2 normalization ==========================================

## Make sample Table ----
sampleTable <- data.frame(
  sampleName = taxa_ids,
  fileName = list.files(
    path = count_path,
    pattern = "\\.counts$",
    full.names = FALSE
  )
)

## Load field metadata (sample) information
nam_design <- read.csv("/workdir/mbb262/nam_hybrid_rnaseq/sample_metadata.csv")

# subset columns
nam_design <- nam_design %>% 
  select("sample_name", "type", "cultivar", "age", "tissue")

# Merge metadata
sampleTable <- merge(x = sampleTable, y = nam_design, 
                     by.x = "sampleName", by.y = "sample_name")

# Remove bad samples
# Based off alignments <50% to b73 on script 8
sampleTable <- sampleTable %>% 
  filter(!(sampleName %in% c("MS21R011", "MS21R014", "MS21R018", 
                             "MS21R025", "MS21R033", "MS21R091", "MS21R098")))
dim(sampleTable)

## Make DESeq data set ----
ddsHTSeq <- DESeq2::DESeqDataSetFromHTSeqCount(
  sampleTable = sampleTable,
  directory = count_path,
  design = ~ 1 # <- no design
)

# turn deseq dataset into datatable object
dtRaw <- SummarizedExperiment::assay(ddsHTSeq) %>%
  as.data.frame() %>%
  data.table::setDT(keep.rownames = "gene_id") %>%
  data.table::melt(id.vars = "gene_id", variable.name = "Taxa") %>%
  data.table::dcast(formula = Taxa ~ gene_id)

# add pseudocount of 1 to count matrix
counts_pseudo <- cbind(dtRaw[,1], dtRaw[,-1] + 1)
counts_pseudo[1:5,1:5]

# remember the names for transposing matrix
n <- counts_pseudo$Taxa

# transpose all but the first column (name)
counts_pseudo <- as.data.frame(t(counts_pseudo[,-1]))
colnames(counts_pseudo) <- n
str(counts_pseudo) # Check the column types

# Gather metadata (coldata) and format
pseudo_coldata <- sampleTable
rownames(pseudo_coldata) <- pseudo_coldata[,1]
pseudo_coldata <- pseudo_coldata[,-c(1,2,3)]

# check if inputs are the correct dimensions
ncol(counts_pseudo)
nrow(pseudo_coldata)

# turn that matrix into a deseq object for vst
dds_pseudo <- DESeqDataSetFromMatrix(countData = counts_pseudo,
                                     colData = pseudo_coldata,
                                     design = ~ 1)

## Variance stabilization transformation ----
ddsVST <- DESeq2::vst(dds_pseudo, blind = FALSE)


# === Save data =====================================================

## dcast arrays as `data.table` objects ----

### Raw count data
dtRaw <- SummarizedExperiment::assay(ddsHTSeq) %>%
  as.data.frame() %>%
  data.table::setDT(keep.rownames = "gene_id") %>%
  data.table::melt(id.vars = "gene_id", variable.name = "Taxa") %>%
  data.table::dcast(formula = Taxa ~ gene_id)

### VST count data
dtVST <- SummarizedExperiment::assay(ddsVST) %>%
  as.data.frame() %>%
  data.table::setDT(keep.rownames = "gene_id") %>%
  data.table::melt(id.vars = "gene_id", variable.name = "Taxa") %>%
  data.table::dcast(formula = Taxa ~ gene_id)


## Write raw data to disk ----

### Raw count data
data.table::fwrite(
  x = dtRaw,
  file = "/workdir/mbb262/nam_hybrid_rnaseq/output/counts/nam_counts_raw.csv"
)

### VST count data
data.table::fwrite(
  x = dtVST,
  file = "/workdir/mbb262/nam_hybrid_rnaseq/output/counts/nam_counts_vst.csv"
)



# # === Filter low counts and taxa ====================================
# # for after i map to correct founder genome
# ## Read in bad quality taxa IDs ----
# 
# ### NOTE: This data was calculated in script `08`
# bad_taxa <- readr::read_lines(file = "/workdir/mbb262/nam_hybrid_rnaseq/output/counts/bad_taxa.txt")
# 
# 
# ## Filter taxa ID ----
# 
# # Minimum filtering at this stage
# min_count <- 10
# number_of_samples <- ddsHTSeq %>%
#     assay() %>%
#     dim() %>%
#     .[2]
# compare <- counts(ddsHTSeq) >= min_count
# 
# ddsHTSeq_filt <- ddsHTSeq[, !(colnames(ddsHTSeq) %in% bad_taxa)]
# ddsHTSeq_filt <- ddsHTSeq_filt[rowSums(compare) > number_of_samples * 0.05, ]
# 
# 
# ## Median of ratios normalization ----
# ddsHTSeq_filt_norm <- DESeq2::estimateSizeFactors(ddsHTSeq_filt)
# 
# 
# 
# ## Variance stabilization transformation ----
# ddsVST_filt <- DESeq2::vst(ddsHTSeq_filt, blind = FALSE)
# 
# 
# ## dcast arrays as `data.table` objects ----
# 
# ### Raw count data
# dtRaw_filt <- SummarizedExperiment::assay(ddsHTSeq_filt) %>%
#     as.data.frame() %>%
#     data.table::setDT(keep.rownames = "gene_id") %>%
#     data.table::melt(id.vars = "gene_id", variable.name = "Taxa") %>%
#     data.table::dcast(formula = Taxa ~ gene_id)
# 
# ### VST count data
# dtVST_filt <- SummarizedExperiment::assay(ddsVST_filt) %>%
#     as.data.frame() %>%
#     data.table::setDT(keep.rownames = "gene_id") %>%
#     data.table::melt(id.vars = "gene_id", variable.name = "Taxa") %>%
#     data.table::dcast(formula = Taxa ~ gene_id)
# 
# ### Normalized count data
# dtNorm_filt <- counts(ddsHTSeq_filt_norm, normalized = TRUE) %>%
#     as.data.frame() %>%
#     data.table::setDT(keep.rownames = "gene_id") %>%
#     data.table::melt(id.vars = "gene_id", variable.name = "Taxa") %>%
#     data.table::dcast(formula = Taxa ~ gene_id)
# 
# 
# ## Write to disk ----
# 
# ### Raw count data
# data.table::fwrite(
#     x = dtRaw_filt,
#     file = "/workdir/mbb262/nam_hybrid_rnaseq/output/counts/nam_counts_raw_filt.csv"
# )
# 
# ### VST count data
# data.table::fwrite(
#     x = dtVST_filt,
#     file = "/workdir/mbb262/nam_hybrid_rnaseq/output/counts/nam_counts_vst_filt.csv"
# )
# 
# ### Normalized count data
# data.table::fwrite(
#     x = dtNorm_filt,
#     file = "/workdir/mbb262/nam_hybrid_rnaseq/output/counts/nam_counts_norm_filt.csv"
# )
# 
# 
