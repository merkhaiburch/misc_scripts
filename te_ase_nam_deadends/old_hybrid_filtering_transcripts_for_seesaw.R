# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-11-27
# Updated... 2023-11-27
#
# Description:
# Old code from hybrids to filter quant.sf files down to the same set of transcripts
# dropping this because the bootstrap files also need to be filtered to the same
# set of bootstrap files but it would require importing and manipulating binary
# files which I wanted to avoid.
# ------------------------------------------------------------------------------

## Find the same set of transcripts across quant files ----------------------------------------------

# tximport relies on the dimensions of the quant.sf files being the same as the 
# transcript/gene file (tx2gene), because of the steps using liftoff, not all B73
# genes will be present in the other founder genomes. This horrifying loop reads
# in all quant.sf files, finds the common set of transcripts found across all
# genomes, reloads the genomes back into r, subsets the quant.sf file to the 
# shared set of transcripts, and exports it as quant_sharedTranscripts.sf

# Make an empty list that will hold all transcript names
allTranscriptIDs <- list()

# Make an empty df that will be used to count the number of transcripts to filter by later
allTranscriptIDs_df <- data.frame()

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
  allTranscriptIDs_df <- rbind(allTranscriptIDs_df, data.frame(temp, id = rep(names(files)[i])))
  rm(temp)
}

# Get a vector of transcripts that are only present in one genome
allTranscriptIDs_df$Name2 <- gsub("(T[0-9]+).*", "\\1", allTranscriptIDs_df$Name)
lala <- allTranscriptIDs_df %>% group_by(Name2) %>% summarise(n = n())

# Remove the source genome
allTranscriptIDs <- lapply(allTranscriptIDs, function(x) gsub("(T[0-9]+).*", "\\1", x))

# Get the common set of transcript ids measured across all samples
tx2gene <- Reduce(base::intersect, allTranscriptIDs) %>% data.frame()
colnames(tx2gene) <- "tx_id"
head(tx2gene)
dim(tx2gene) # 69,238 

# Turn into a df, column 1 transcripts, column 2 gene name
# tx2gene is a three-column dataframe linking transcript ID (column 1) to 
# gene ID (column 2). We will take the first two columns as input to tximport. 
tx2gene$geneID <- gsub("_.*", "", tx2gene$tx_id)
tx2gene <- tx2gene[order(tx2gene$tx_id),] # sort alphabetically
head(tx2gene)

# Need to also intersect the above ids with what's in the gff file
# Get coordinates for transcripts from gff
# wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz
b73_gff <- ape::read.gff("/workdir/mbb262/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3") %>% 
  dplyr::filter(type == "mRNA" & seqid %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10")) %>% 
  dplyr::select("seqid", "start", "end", "strand", "attributes")
b73_gff$tx_id <- gsub("ID=", "", b73_gff$attributes)
b73_gff$tx_id <- gsub(";.*", "", b73_gff$tx_id)
b73_gff$geneID <- gsub(".*Parent=", "", b73_gff$attributes)
b73_gff$geneID <- gsub(";.*", "", b73_gff$geneID)
b73_gff <- b73_gff %>% dplyr::select(-"attributes")
head(b73_gff)

# Filter to only those transcripts measured across all samples AND that are present in the gff
temp <- dplyr::inner_join(tx2gene, b73_gff, by = "tx_id") %>% select(-"geneID.y")
colnames(temp)[2] <- "geneID"
tx2gene <- temp %>% dplyr::select(tx_id, geneID)
b73_gff <- temp %>% dplyr::select(seqid, start, end, strand, tx_id, geneID)
colnames(b73_gff)[6] <- "gene_id" # chamge column name into what fishpond will accept

# Check dimensions, should be the same, 68,969
dim(temp)
dim(b73_gff)
dim(tx2gene)
head(b73_gff)
head(tx2gene)

# Turn into a genomic ranges object
combined_gff_gr <- GenomicRanges::makeGRangesFromDataFrame(b73_gff, keep.extra.columns = TRUE)

# Link transcripts with genes (tx2gene)
t2g <- fishpond::makeTx2Tss(combined_gff_gr) # GRanges object
GenomicRanges::mcols(t2g)[,c("tx_id","group_id")]


## Filter salmon files to the same consistent set of transcripts ------------------------------------------

# Load more problematic packages
library(GenomicRanges)
library(fishpond)

# Loop through all files, subset to the shared 68,398 transcripts, export with a different name
for (i in 1:length(samples)){
  # Read in file
  temp <- data.table::fread(files[i], nThread = 40)
  
  # Make new column that removes the genome source name
  temp$shortTranscript <- gsub("(T[0-9]+).*", "\\1", temp$Name)
  
  # # Filter out transcripts to only those shared across all files,
  # # then remove that column
  # temp <- temp %>%
  #   dplyr::filter(shortTranscript %in% tx2gene$tx_id) %>%
  #   dplyr::select(-shortTranscript)
  
  # Make new column that removes the genome source name
  temp$shortTranscript <- gsub("(T[0-9]+).*", "\\1", temp$Name)
  
  # Filter out transcripts that are seen more than once
  lala <- temp %>%
    group_by(shortTranscript) %>% 
    summarise(n = n()) %>% 
    dplyr::filter(n > 1)
  temp <- temp %>% 
    dplyr::filter(shortTranscript %in% lala$shortTranscript) %>% 
    dplyr::select(-shortTranscript)
  
  # Check dimensions
  # For hybrids, its normal that these counts aren't consistent between dfs,
  # some count a transcript more than once
  print(dim(temp))
  print(head(temp))
  
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
