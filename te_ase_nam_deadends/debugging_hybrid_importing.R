# Debugging hybrid imports

# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-11-08
# Updated... 2023-11-08
#
# Description:
# Get hybrid expression and allele specific expression using SeeSaw
# https://mikelove.github.io/Bioc2022AllelicExpression/articles/Bioc2022AllelicExpression.html
# ------------------------------------------------------------------------------

# Load packages
library(dplyr)


## Key error where this all begins ---------------------------------------------

# Find files
samples <- list.files("/workdir/mbb262/nam_hybrid_rnaseq/salmon/output/hybrid", full.names = TRUE)
files <- file.path(samples, "quant.sf")
# Subset to debug --> c("MS21R193", "MS21R254", "MS21R317", "MS21R384")
samples <- samples[c(60,120, 183, 248)]
files <- files[c(60,120, 183, 248)]
shared_files <- file.path(samples, "shared_quant.sf")
names(shared_files) <- gsub(".*/", "", samples)
shared_files_df <- data.frame(shared_files)
shared_files_df$sample_name <- rownames(shared_files_df)

# Make metadata for loading into R
meta_hybrids <- read.csv("/home/mbb262/git_projects/te_ase_nam/data/sample_metadata.csv") %>% 
  dplyr::filter(type == "hybrid") %>% 
  dplyr::select(sample_name, cultivar, group, age, tissue)
temp <- table(meta_hybrids$cultivar)
coldata <- merge(meta_hybrids, shared_files_df, by = "sample_name") 
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
                         a1="alt", 
                         a2="ref",
                         format="assays", 
                         tx2gene=t2g,
                         dropInfReps = FALSE)

# Error message
# reading in files with read_tsv
# 1 2 3 4 
# transcripts missing from tx2gene: 414 --> successfully recreated this
# summarizing abundance --> successfully recreated this
# summarizing counts --> successfully recreated this
# summarizing length --> successfully recreated this
# Error in rowsum.default(x[sub.idx, , drop = FALSE], geneId) : 
#   incorrect length for 'group'

## Debugging hybrid read import ------------------------------------------------

# from tximport 
varReduce <- FALSE #default
txOut <- FALSE #default
infRepType <- if (varReduce & txOut) { "var" } else { "full" }

# Key piece of code for other variables
txIdCol <- "Name"
abundanceCol <- "TPM"
countsCol <- "NumReads"
lengthCol <- "EffectiveLength"

for (i in seq_along(shared_files)) {
  message(i)
  
  # Load in inf reps
  repInfo <- readInfRepFish(dirname(shared_files[i]))
  
  # import and convert quantification info to data.frame
  raw <- read.delim(shared_files[i])
  
  if (i == 1) {
    txId <- raw[[txIdCol]]
  } else {
    stopifnot(all(txId == raw[[txIdCol]]))
  }
  
  if (i == 1) {
    mat <- matrix(nrow=nrow(raw),ncol=length(files))
    rownames(mat) <- raw[[txIdCol]]
    colnames(mat) <- names(files)
    abundanceMatTx <- mat
    countsMatTx <- mat
    lengthMatTx <- mat
    
    varMatTx <- mat
    
    infRepMatTx <- list()
  }
  
  abundanceMatTx[,i] <- raw[[abundanceCol]]
  countsMatTx[,i] <- raw[[countsCol]]
  lengthMatTx[,i] <- raw[[lengthCol]]
  infRepMatTx[[i]] <- repInfo$reps
  
  # repInfo <- tximport::.infRepImporter(dirname(files[1]))
  # varMatTx[,i] <- repInfo$vars
}


txi <- list(abundance=abundanceMatTx,
            counts=countsMatTx,
            length=lengthMatTx,
            infReps=infRepMatTx)

# txi[["countsFromAbundance"]] <- NULL
# txiGene <- summarizeToGene(txi, tx2gene, varReduce, ignoreTxVersion, ignoreAfterBar, countsFromAbundance)

# unpack matrices from list for cleaner code --> summarizeToGene code from tximport
# https://rdrr.io/bioc/tximport/src/R/summarizeToGene.R
abundanceMatTx <- txi$abundance
countsMatTx <- txi$counts
lengthMatTx <- txi$length
txId <- rownames(abundanceMatTx)
txId <- sub("_alt", "", txId) # remove everything after underscore (_ref, _alt)
txId <- sub("_ref", "", txId) # remove everything after underscore (_ref, _alt)

# I added this to simplify debugging
tx2gene <- t2g

# Code from helper_allelic.R
cols <- c("tx_id","group_id")
# swap around variable names to run the data.frame-based code
txps <- tx2gene
tx2gene <- GenomicRanges::mcols(txps)[,cols]

# back to summarizeToGene code from tximport
# https://rdrr.io/bioc/tximport/src/R/summarizeToGene.R
colnames(tx2gene) <- c("tx","gene")
tx2gene$gene <- factor(tx2gene$gene)
tx2gene$tx <- factor(tx2gene$tx)
tx2gene

# remove transcripts (and genes) not in the rownames of matrices
tx2gene <- tx2gene[tx2gene$tx %in% txId,]
tx2gene$gene <- droplevels(tx2gene$gene)
ntxmissing <- sum(!txId %in% tx2gene$tx)
if (ntxmissing > 0) message("transcripts missing from tx2gene: ", ntxmissing) # matches up with printed error message abouve

# subset to transcripts in the tx2gene table
sub.idx <- txId %in% tx2gene$tx
abundanceMatTx <- abundanceMatTx[sub.idx,,drop=FALSE]
countsMatTx <- countsMatTx[sub.idx,,drop=FALSE]
lengthMatTx <- lengthMatTx[sub.idx,,drop=FALSE]
txId <- txId[sub.idx]

# now create a vector of geneId which aligns to the matrices
geneId <- tx2gene$gene[match(txId, tx2gene$tx)]

# summarize abundance and counts
message("summarizing abundance")
abundanceMat <- rowsum(abundanceMatTx, geneId)
message("summarizing counts")
countsMat <- rowsum(countsMatTx, geneId)
message("summarizing length")

# Run the key element that's been messing with me
if ("infReps" %in% names(txi)) {
  infReps <- lapply(txi$infReps, function(x) rowsum(x[sub.idx,,drop=FALSE], geneId))
  message("summarizing inferential replicates")
}

# Try running on just one element of the list
test <- rowsum(x[sub.idx,,drop=FALSE], geneId)
x <- txi$infReps[[1]]
xsub <- x[sub.idx,,drop=FALSE]
dim(xsub) #139821
length(geneId) #137014


# Traceback call
# Error in rowsum.default(x[sub.idx, , drop = FALSE], geneId) :
#   incorrect length for 'group'
# 5. stop("incorrect length for 'group'")
# 4. rowsum.default(x[sub.idx, , drop = FALSE], geneId)
# 3. rowsum(x[sub.idx, , drop = FALSE], geneId)
# 2. FUN(X[[i]], ...)
# 1. lapply(txi$infReps, function(x) rowsum(x[sub.idx, , drop = FALSE], geneId))

# from debug
function (x, group, reorder = TRUE, na.rm = FALSE, ...) 
{
  if (!is.numeric(x)) 
    stop("'x' must be numeric")
  if (length(group) != NROW(x)) 
    stop("incorrect length for 'group'") # where its messing up
  if (anyNA(group)) 
    warning("missing values for 'group'")
  ugroup <- unique(group)
  if (reorder) 
    ugroup <- sort(ugroup, na.last = TRUE, method = "quick")
  .Internal(rowsum_matrix(x, group, ugroup, na.rm, as.character(ugroup)))
}





# End debugging hybrid read import ---------------------------------------------
































## Find files ------------------------------------------------------------------

# List all directories containing data 
samples <- list.files("/workdir/mbb262/nam_hybrid_rnaseq/salmon/output/hybrid", full.names = TRUE)

# Obtain a vector of all filenames including the path
files <- file.path(samples, "quant.sf")

# Since all quant files have the same name it is useful to have names for each element
names(files) <- stringr::str_replace(samples, ".*/hybrid/", "") 

# Subset to just three files to debug this pipeline
files <- files[1:3]
samples <- samples[1:3]


## Find the same set of transcripts across quant files ----------------------------------------------

# Make an empty list that will hold all transcript names
allTranscriptIDs <- list()

# Run loop
for (i in 1:length(samples)){
  # Read in file
  temp <- data.table::fread(files[i], nThread = 40) %>% 
    dplyr::select("Name") %>% 
    as.list()
  
  # Print the sample I'm on
  print(names(files[i]))
  
  # Get names of transcripts
  allTranscriptIDs <- c(allTranscriptIDs, temp)
  rm(temp)
}

# Remove the source genome
allTranscriptIDs <- lapply(allTranscriptIDs, function(x) gsub("(T[0-9]+).*", "\\1", x))

# Get the common set of transcript ids measured across all samples
tx2gene <- Reduce(base::intersect, allTranscriptIDs) 
length(tx2gene) # 70534

# Turn into a df, column 1 transcripts, column 2 gene name
tx2gene <- data.frame(tx2gene)
tx2gene$geneID <- gsub("_.*", "", tx2gene$tx2gene)
colnames(tx2gene)[1] <- "Name"
tx2gene <- tx2gene[order(tx2gene$Name),] # sort alphabetically
head(tx2gene)
temp <- table(tx2gene$geneID)


# Need to also intersect the above ids with what's in the gff file
# Get coordinates for transcripts from gff
# wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz
b73_gff <- ape::read.gff("/workdir/mbb262/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3") %>% 
  dplyr::filter(type == "mRNA" & seqid %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10")) %>% 
  dplyr::select("seqid", "start", "end", "strand", "attributes")
b73_gff$tx_id <- gsub("ID=", "", b73_gff$attributes)
b73_gff$tx_id <- gsub(";.*", "", b73_gff$tx_id)
b73_gff$gene_id <- gsub(".*Parent=", "", b73_gff$attributes)
b73_gff$gene_id <- gsub(";.*", "", b73_gff$gene_id)
b73_gff <- b73_gff %>% dplyr::select(-"attributes")
colnames(b73_gff)[5] <- "Name"

# Filter to only those measured across all samples AND that are present in the gff
temp <- dplyr::inner_join(tx2gene, b73_gff, by = "Name")
tx2gene <- temp %>% dplyr::select(Name, geneID)
b73_gff <- temp %>% dplyr::select(seqid, start, end, strand, Name, gene_id)
colnames(b73_gff)[5] <- "tx_id" # change column name back

# Check dimensions, should be the same --> 70,266
dim(b73_gff)
dim(tx2gene)

# Turn into a genomic ranges object
combined_gff_gr <- GenomicRanges::makeGRangesFromDataFrame(b73_gff, keep.extra.columns = TRUE)

# Link transcripts with genes (tx2gene)
t2g <- fishpond::makeTx2Tss(combined_gff_gr) # GRanges object
GenomicRanges::mcols(t2g)[,c("tx_id","group_id")]


## Filter salmon files to the same consistent set of transcripts ------------------------------------------

# Loop through all files, subset to the shared 68,398 transcripts, export with a different name
for (i in 1:length(samples)){
  # Read in file
  temp <- data.table::fread(files[i], nThread = 40)
  
  # Make new column that removes the genome source name
  temp$shortTranscript <- gsub("(T[0-9]+).*", "\\1", temp$Name)
  
  # Filter out transcripts to only those shared across all files,
  # then remove that column
  temp <- temp %>%
    dplyr::filter(shortTranscript %in% tx2gene$Name) %>%
    dplyr::select(-shortTranscript)
  
  # Check dimensions
  # For hybrids, its normal that these counts aren't consistent between dfs,
  # some count a transcript more than once
  print(dim(temp))
  
  # Sort alphabetically because of this error
  # https://support.bioconductor.org/p/96588/
  temp <- temp[order(temp$Name),]
  
  # Simplify names for fishpond to ref or alt
  temp$Name <- sub("B73|PHB47|PHZ51", "ref", temp$Name)
  temp$Name <- sub("B97|CML52|CML69|CML103|CML228|CML247|CML277|CML322|CML333|HP301|Il14H|Ki3|Ki11|Ky21|M37W|M162W|Mo17|Mo18W|Ms71|NC350|NC358|Oh7B|Oh43|P39|Tx303|Tzi8", 
                   "alt", temp$Name)
  
  # Save with a different name
  updated_name <- gsub("quant.sf", "shared_quant.sf", files)
  data.table::fwrite(temp, file = updated_name[i], nThread = 40, sep = "\t")
  rm(temp)
  gc()
}

# Obtain a vector of all modified filenames including the path
shared_files <- file.path(samples, "shared_quant.sf")
names(shared_files) <- gsub(".*/", "", samples)

###################################################################################

# Debug room

# load files
file1 <- data.table::fread(files[1], nThread = 40)
file2 <- data.table::fread(files[2], nThread = 40)
file3 <- data.table::fread(files[3], nThread = 40)

# Remove genome source
file1$txid <- gsub("(T[0-9]+).*", "\\1", file1$Name)
file2$txid <- gsub("(T[0-9]+).*", "\\1", file2$Name)
file3$txid <- gsub("(T[0-9]+).*", "\\1", file3$Name)

# Count the number of unique transcripts 
length(unique(tx2gene$Name)) # 70266
dim(tx2gene) # also the same

# Remove duplicates
file1 <- file1 %>% dplyr::filter(txid %in% tx2gene$Name)
file2 <- file2 %>% dplyr::filter(txid %in% tx2gene$Name)
file3 <- file3 %>% dplyr::filter(txid %in% tx2gene$Name)

# Count dimensions --> they're all different and they shouldn't be, what's doubled up?
dim(file1)
dim(file2)
dim(file3)
length(unique(file1$txid)) #70266 as it should be....
length(unique(file2$txid)) #70266 as it should be....
length(unique(file3$txid)) #70266 as it should be....

# Count the number of transcripts, each transcript should = 2
temp1 <- file1 %>% group_by(txid) %>% summarize(n = n())
temp2 <- file2 %>% group_by(txid) %>% summarize(n = n())
temp3 <- file3 %>% group_by(txid) %>% summarize(n = n())

# Count how many transcripts = 1 (should be 2 if its shared in both genomes)
length(which(temp1$n != 1)) # 60021
length(which(temp2$n != 1)) # 60985
length(which(temp3$n != 1)) # 59755

# Try combining all and recounting
lala <- rbind(file1, file2, file3)
temp4 <- lala %>% group_by(txid) %>% summarize(n = n())

#test intersect with list
lala <- Reduce(base::intersect,  list(v1 = file1$txid, 
                        v2 = file2$txid, 
                        v3 = file3$txid))
length(lala)


###################################################################################
## Create SEESAW friendly metadata ---------------------------------------------

# Load more problematic packages
library(GenomicRanges)
library(fishpond)

# Turn file paths to a df
shared_files_df <- data.frame(shared_files)
shared_files_df$sample_name <- rownames(shared_files_df)
dim(shared_files_df)

# Load metadata
meta_hybrids <- read.csv("/home/mbb262/git_projects/te_ase_nam/data/sample_metadata.csv") %>% 
  dplyr::filter(type == "hybrid") %>% 
  dplyr::select(sample_name, cultivar, group, age, tissue)

# Merge metadata with paths
coldata <- merge(meta_hybrids, shared_files_df, by = "sample_name")
dim(coldata)
head(coldata)


## Load in salmon quants using fishpond ----------------------------------------

# Split importing reads into three groups, b73, phz51, and phb47
gse_b73 <- importAllelicCounts(coldata, 
                               a1="alt", 
                               a2="ref",
                               format="wide", 
                               tx2gene=t2g,
)


# `files` points to `quant.sf` files in original Salmon directories
files <- fishpond::importAllelicCounts("/workdir/mbb262/sim_reads/salmon/hybridb73ky21/hybrid_transcripts_quant/shared_quant.sf",
                                       a1 = 'ky21', a2 = 'b73', tx2gene = t2g, dropInfReps = TRUE)

SummarizedExperiment::colData(files)
files <- labelKeep(files)
files <- swish(files, x="allele", pair="sample", cor="pearson")


## More debugging stuff from another script ------------------------------------

##################################################################################
# DEBUGGING SECTION

temp <- read.delim(shared_files[1])
temp$shortTranscript <- gsub("(T[0-9]+).*", "\\1", temp$Name)

lala <- temp %>%
  group_by(shortTranscript) %>% 
  summarise(n = n()) %>% 
  dplyr::filter(n > 1)

temp <- temp %>% 
  dplyr::filter(shortTranscript %in% lala$shortTranscript) %>% 
  dplyr::select(-shortTranscript)

rownames(temp) <- temp$Name
temp$shortName <- gsub("_alt", "", temp$Name)
temp$shortName <- gsub("_ref", "", temp$shortName)

# Debug some shit
(ntxp <- nrow(temp)/2)
(n <- ncol(temp))

a1match <- paste0("_", "alt","$")
a2match <- paste0("_", "ref","$")

# gather transcript names for a1 and a2 alleles
txp_nms_a1 <- grep(a1match, rownames(temp), value=TRUE)
length(txp_nms_a1) == ntxp

# from tximport ----------------------------------------------

# Key piece of code for other variables
txIdCol <- "Name"
abundanceCol <- "TPM"
countsCol <- "NumReads"
lengthCol <- "EffectiveLength"

# abundanceMatTx[,i] <- raw[[abundanceCol]]
# countsMatTx[,i] <- raw[[countsCol]]
# lengthMatTx[,i] <- raw[[lengthCol]]

i <- 1
for (i in seq_along(shared_files)) {
  
  # import and convert quantification info to data.frame
  raw <- read.delim(shared_files[i])
  
  if (i == 1) {
    txId <- raw[[txIdCol]]
  } else {
    stopifnot(all(txId == raw[[txIdCol]]))
  }
  
  if (i == 1) {
    mat <- matrix(nrow=nrow(raw),ncol=length(files))
    rownames(mat) <- raw[[txIdCol]]
    colnames(mat) <- names(files)
    abundanceMatTx <- mat
    countsMatTx <- mat
    lengthMatTx <- mat
    
    varMatTx <- mat
    # if (infRepType == "var") {
    #   varMatTx <- mat
    # } else if (infRepType == "full") {
    #   infRepMatTx <- list()
    # }
  }
  abundanceMatTx[,i] <- raw[[abundanceCol]]
  countsMatTx[,i] <- raw[[countsCol]]
  lengthMatTx[,i] <- raw[[lengthCol]]
  
  # repInfo <- tximport::.infRepImporter(dirname(files[1]))
  # varMatTx[,i] <- repInfo$vars
}


txi <- list(abundance=abundanceMatTx,
            counts=countsMatTx,
            length=lengthMatTx)

# txi[["countsFromAbundance"]] <- NULL
# txiGene <- summarizeToGene(txi, tx2gene, varReduce, ignoreTxVersion, ignoreAfterBar, countsFromAbundance)

# unpack matrices from list for cleaner code --> summarizeToGene code
abundanceMatTx <- txi$abundance
countsMatTx <- txi$counts
lengthMatTx <- txi$length
txId <- rownames(abundanceMatTx)

# Mer added
tx2gene <- t2g # mer added

# Code from helper_allelic.R
cols <- c("tx_id","group_id")
# swap around variable names to run the data.frame-based code
txps <- tx2gene
tx2gene <- mcols(txps)[,cols]


# back to summarizeToGene code
colnames(tx2gene) <- c("tx","gene")
tx2gene$gene <- factor(tx2gene$gene)
tx2gene$tx <- factor(tx2gene$tx)
tx2gene

# remove transcripts (and genes) not in the rownames of matrices
tx2gene <- tx2gene[tx2gene$tx %in% txId,]
tx2gene$gene <- droplevels(tx2gene$gene)
ntxmissing <- sum(!txId %in% tx2gene$tx)
if (ntxmissing > 0) message("transcripts missing from tx2gene: ", ntxmissing)

# subset to transcripts in the tx2gene table
sub.idx <- txId %in% tx2gene$tx
abundanceMatTx <- abundanceMatTx[sub.idx,,drop=FALSE]
countsMatTx <- countsMatTx[sub.idx,,drop=FALSE]
lengthMatTx <- lengthMatTx[sub.idx,,drop=FALSE]
txId <- txId[sub.idx]

# now create a vector of geneId which aligns to the matrices
geneId <- tx2gene$gene[match(txId, tx2gene$tx)]


# Old code -----------------------------------------------------------------------------------

# Load data
hybrid_vst <- read.delim()

# Load sample metadata
sample_metadata <- read.csv("/home/mbb262/git_projects/te_ase_nam/")

# Old code, haven't changed from test example yet

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

# Old code -----------------------------------------------------------------------------------

# Load data
hybrid_vst <- read.delim()

# Load sample metadata
sample_metadata <- read.csv("/home/mbb262/git_projects/te_ase_nam/")

# Old code, haven't changed from test example yet

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

