# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-09-13 
# Updated... 2022-09-15

# Description 
# Run a loop over each gene and calculate kinship matrix, apply it
# and fit a mixed model
# expression_i ~ TE PAV_i + K_maize (where i = each individual gene)
# ------------------------------------------------------------------------------

# Restart R twice

# Set memory options
options(java.parameters = "-Xmx150g")
options(java.parameters = "-Xmx150g")
options(java.parameters = "-Xmx150g")
library(rJava)
.jinit()
memUsageAfter <- rJava::.jcall(rJava::.jnew("java/lang/Runtime"), "J", "maxMemory") / 1024^3
message("Your max memory usage after changing it is: ", memUsageAfter, "GB")

# Load in packages
library(rTASSEL)
library(dplyr)
library(data.table)

# Start logger
rTASSEL::startLogger(fullPath = "/workdir/mbb262", fileName = "rtassel_logger.txt")

# Set the number of threads
numThreads <- 50


# Load in maize vcf file to calculate K with -----------------------------------
# Load in genotype table for inbreds and hybrids
# of SNPs: 10,156,186; of taxa: 113
# mkdir /workdir/mbb262/genotypes
# scp -r mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/mbb262/te_ase_nam/genotypes_for_k/v5_inbred_hybrid_genotypes/hmp321_282_agpv5_inbred_hybrid_all_chroms.vcf.gz /workdir/mbb262/genotypes
tasGeno <-  rTASSEL::readGenotypeTableFromPath("/workdir/mbb262/genotypes/hmp321_282_agpv5_inbred_hybrid_all_chroms.vcf.gz", keepDepth = FALSE)
print(tasGeno)



# Format TE PAV matrix ---------------------------------------------------------

# Download file from CBSU
# scp mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/mcs368/anchorwave_indel/NAM_anchorwave_indels.PAV01.noAncestralSpecification.2022-09-27.txt.gz /workdir/mbb262/genotypes
# gunzip /workdir/mbb262/genotypes/NAM_anchorwave_indels.PAV01.noAncestralSpecification.2022-09-27.txt.gz

# Load file into R
temp_te <- data.table::fread("/workdir/mbb262/genotypes/NAM_anchorwave_indels.PAV01.noAncestralSpecification.2022-09-27.txt", 
                             sep = " ")

# transpose matrix
temp_te <- t(temp_te) %>% data.frame()

# Create marker names
colnames(temp_te) <- paste0("te_", stringr::str_pad(1:ncol(temp_te), 7, pad = "0"))

# Save first two rows to cross-reference later, remove them from TE file
te_key_file <- t(temp_te[1:2,]) %>% data.frame()
te_key_file$markerID <- rownames(te_key_file)
temp_te <- temp_te[-c(1:2),]

# Change rownames to hapmap names given a handmade cross-reference file
x_ref_nam_hapmap <- read.csv("/home/mbb262/git_projects/te_ase_nam/data/nam_hapmap_xref.csv")

# Separeate out rownames
temp_te$rowID <- rownames(temp_te)

# Merge x-ref file with table
temp_te <- merge(temp_te, x_ref_nam_hapmap, by.x = "rowID", by.y = "NAM")
rownames(temp_te) <- temp_te$hapmap321
temp_te <- temp_te %>% select(-"rowID", -"hapmap321") # remove these extra column now

# Check
temp_te[1:3,1:3]

# Write to file
data.table::fwrite(temp_te, "/workdir/mbb262/genotypes/NAM_anchorwave_indels.PAV01.noAncestralSpecification.2022-09-27_rformat.txt",
                   sep = "\t", row.names = TRUE)

# Quick and dirty trick, open in sublime or vim, add these two lines to top of document
# <Numeric>
#   <Marker>

#############################
# Use the following script to create hybrid genotypes: model_create_hybrid_te_genotypes.sh
#############################

# Load in FORMATTED TE PAV matrix for inbreds and hybrids
te_genotypes_inbreds <-  data.table::fread("/workdir/mbb262/genotypes/hybrid_numeric_geno_no_header.txt", nThread = numThreads)
colnames(te_genotypes_inbreds)[1] <- "Taxa"

# Just subset to inbred parents for test
# te_genotypes_inbreds <- temp_te[1:26,]


# Load in expression -----------------------------------------------------------

# Load in expression data (taxa x genes)

# For now simulate it
expression_df <- x_ref_nam_hapmap$hapmap321 %>% data.frame()
colnames(expression_df) <- "Taxa"
expression_df$Zm00001eb000050 <- rnorm(n = nrow(expression_df))
expression_df$Zm00001eb000080 <- rnorm(nrow(expression_df))
expression_df$Zm00001eb000170 <- rnorm(nrow(expression_df))
expression_df$Zm00001eb000180 <- rnorm(nrow(expression_df))
expression_df$Zm00001eb000190 <- rnorm(nrow(expression_df))
expression_df$Zm00001eb000210 <- rnorm(nrow(expression_df))
expression_df$Zm00001eb000210 <- rnorm(nrow(expression_df))
expression_df$Zm00001eb000240 <- rnorm(nrow(expression_df))
expression_df$Zm00001eb000310 <- rnorm(nrow(expression_df))
expression_df$Zm00001eb000330 <- rnorm(nrow(expression_df))



# Create windows around genes for maize K matrices -----------------------------

# Download files
# system("wget -P /workdir/mbb262/ https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz")
# system("wget -P /workdir/mbb262/ https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.canonical_transcripts")
# system("wget -P /workdir/mbb262/ https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.canonical_transcripts")

# Load in gff file, filter to 5' UTR regions
gff_5utr <- ape::read.gff("/workdir/mbb262/k_files/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz", na.strings = c(".", "?"), GFF3 = TRUE) %>% 
  filter(type == "five_prime_UTR", seqid %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10")) %>% 
  select("seqid", "start", "end", "strand", "attributes")

# Parse out gene names from metadata column
gff_5utr$attributes <- gsub("Parent=", "", gff_5utr$attributes)

# Subset gff just canonical transcripts
canonical <- read.delim("/workdir/mbb262/k_files/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.canonical_transcripts", header = FALSE)
gff_5utr <- gff_5utr %>% filter(attributes %in% canonical$V1) 

# Gather the transcript with the smallest start position
gff_5utr <- gff_5utr %>% 
  group_by(attributes) %>% 
  slice(which.min(start))

# add a fixed window size to each interval
window_size <- 500
gff_5utr$start <- gff_5utr$start - window_size

# Remove _TXXX from genes
gff_5utr$attributes <- gsub("_T[0-9]{3}", "", gff_5utr$attributes)



# Annotate TE genotype file with genes they're nearest -------------------------

# Label what TEs are closest to genes

# Get and format gene gff file
gff_gene <- ape::read.gff("/workdir/mbb262/k_files/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz", 
                          na.strings = c(".", "?"), GFF3 = TRUE) %>% 
  filter(type == "gene", seqid %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10")) %>% 
  select("seqid", "start", "end", "strand", "attributes")
gff_gene$attributes <- gsub("ID=", "", gff_gene$attributes) # parse out gene names
gff_gene$attributes <- gsub(";.*", "", gff_gene$attributes) # parse out gene names

# Window to call a TE close to a gene
te_window <- 5000

# Add a window to gene's start
gff_gene$start <- gff_gene$start - te_window

# Format the TE marker file
te_key_file$CHROM <- paste0("chr", te_key_file$CHROM)
rownames(te_key_file) <- NULL
te_key_file$end <- te_key_file$POS
colnames(te_key_file) <- c("seqid", "start", "markerID", "end")

# Find overlapping TEs within a "te_window" size of a gene
gff_gene <- data.table::as.data.table(gff_gene)
te_key_file <- data.table::as.data.table(te_key_file)
setkey(gff_gene, seqid, start, end)
gene_te_overlaps <- data.table::foverlaps(te_key_file, gff_gene, type="within",nomatch=0L) #0L drops TEs with no overlapping genes

# Select just the TE columns
gene_te_overlaps <- gene_te_overlaps %>% select("attributes", "markerID")

# Fix column names
colnames(gene_te_overlaps) <- c("v5_gene", "te_marker_id")



# Run mixed models -------------------------------------------------------------

# Get list of genes to iterate through that:
# Are R2<0.9 correlated with TE K matrices 
# should save somewhere on blfs1
test_genes <- read.delim("/workdir/mbb262/k_files/kinship_correlation_500bpgene_5kbte.txt")
test_genes <- test_genes %>% filter(Pearson.Correlation^2 <= 0.9) %>% select(Gene.Name)
test_genes <- test_genes$Gene.Name

# Create place to store results
te_model_effects_results <- c()
te_model_stats_results <- c()

# Run loop
start_time_foreach <- Sys.time()

# length(test_genes)
for (gene in 1:5) {
  
  # Subset out single gene's coordinates
  maize_coord <- gff_5utr %>% filter(attributes == test_genes[gene])
  print(test_genes[gene])
  
  # Subset vcf file
  gr_maize_k <- GenomicRanges::makeGRangesFromDataFrame(maize_coord)
  maize_k_vcf <- tasGeno %>% rTASSEL::filterGenotypeTableSites(siteRangeFilterType = "none", gRangesObj = gr_maize_k)
  
  # Calculate maize K
  k_maize <- kinshipMatrix(tasObj = maize_k_vcf)

  # Filter TE genotypes to those near the [gene] of interest
  gene_tes <- gene_te_overlaps %>% data.frame()%>% filter(v5_gene == test_genes[gene]) %>% select(te_marker_id)
  filter_te_geno <- te_genotypes_inbreds[, gene_tes[,1]]
  filter_te_geno <- tibble::rownames_to_column(filter_te_geno, "<Marker>")
  data.table::fwrite(filter_te_geno, "/workdir/mbb262/temp_genotype.txt", sep = "\t")
  system("echo '<Numeric>' | cat - /workdir/mbb262/temp_genotype.txt > temp && mv temp /workdir/mbb262/temp_genotype.txt")
  # system("sed -i '1s/^/<Numeric>\n/'  /workdir/mbb262/temp_genotype.txt)
  # exe = paste0("sed -i '1s/^/<Numeric>\n", "/'")

  # Filter expression matrix to just that gene
  filtered_expression_df <- expression_df %>% select(Taxa, test_genes[gene])
  
  # Turn genotypes and phenotypes into a tassel object
  tasGenoPheno <- rTASSEL::readGenotypePhenotype(
    genoPathOrObj = "/workdir/mbb262/temp_genotype.txt",
    phenoPathDFOrObj = filtered_expression_df,
    taxaID = "Taxa",
    attributeTypes = NULL
  )
  tasGenoPheno
  
  # Run mixed model
  tasMLM <- rTASSEL::assocModelFitter(
    tasObj = tasGenoPheno,
    formula = . ~ .,
    fitMarkers = TRUE,
    kinship = k_maize,
    fastAssociation = FALSE, 
    maxThreads = numThreads)
  
  # Add a label to mixed model results for this gene
  tasMLM$MLM_Effects$gene <- rep(test_genes[gene], nrow(tasMLM$MLM_Effects))
  tasMLM$MLM_Stats$gene <- rep(test_genes[gene], nrow(tasMLM$MLM_Stats))
  
  # Save results outside of loop
  te_model_effects_results <- rbind(te_model_effects_results, tasMLM$MLM_Effects)
  te_model_stats_results <- rbind(te_model_stats_results, tasMLM$MLM_Stats)
}

# Estimate timing
end_time_foreach <- Sys.time()
(end_time_foreach - start_time_foreach)

# Write to file
data.table::fwrite(mlm_effect_results, "/workdir/mbb262/te_model_effects_results", nThread = numThreads)
data.table::fwrite(mlm_stats_results, "/workdir/mbb262/te_model_stats_results", nThread = numThreads)

# Plotting ---------------------------------------------------------------------

