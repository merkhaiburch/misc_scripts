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
options(java.parameters = "-Xmx500g")
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

# # Only run this section if no formatted matrix exists
# # Download file from CBSU
# # scp mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/mcs368/anchorwave_indel/NAM_anchorwave_indels.PAV01.noAncestralSpecification.2022-09-27.txt.gz /workdir/mbb262/genotypes
# # gunzip /workdir/mbb262/genotypes/NAM_anchorwave_indels.PAV01.noAncestralSpecification.2022-09-27.txt.gz
# 
# # Load file into R
# temp_te <- data.table::fread("/workdir/mbb262/genotypes/NAM_anchorwave_indels.PAV01.noAncestralSpecification.2022-09-27.txt",
#                              sep = " ")
# 
# # transpose matrix
# temp_te <- t(temp_te) %>% data.frame()
# 
# # Create marker names
# colnames(temp_te) <- paste0("te_", stringr::str_pad(1:ncol(temp_te), 7, pad = "0"))
# 
# # Save first two rows to cross-reference later, remove them from TE file
# te_key_file <- t(temp_te[1:2,]) %>% data.frame()
# te_key_file$markerID <- rownames(te_key_file)
# write.csv(te_key_file, "/workdir/mbb262/genotypes/te_key_file.csv", row.names = F, quote = F)
# temp_te <- temp_te[-c(1:2),]
# 
# # Change rownames to hapmap names given a handmade cross-reference file
# x_ref_nam_hapmap <- read.csv("/home/mbb262/git_projects/te_ase_nam/data/nam_hapmap_xref.csv")
# 
# # Separeate out rownames
# temp_te$rowID <- rownames(temp_te)
# 
# # Merge x-ref file with table
# temp_te <- merge(temp_te, x_ref_nam_hapmap, by.x = "rowID", by.y = "NAM")
# rownames(temp_te) <- temp_te$hapmap321
# temp_te <- temp_te %>% select(-"rowID", -"hapmap321") # remove these extra column now
# 
# # Check
# temp_te[1:3,1:3]
# 
# # Write to file
# data.table::fwrite(temp_te, "/workdir/mbb262/genotypes/NAM_anchorwave_indels.PAV01.noAncestralSpecification.2022-09-27_rformat.txt",
#                    sep = "\t", row.names = TRUE)
# 
# # Quick and dirty trick, open in sublime or vim, add these two lines to top of document
# # <Numeric>
# #   <Marker>

#############################
# Use the following script to create hybrid genotypes: model_create_hybrid_te_genotypes.sh
#############################

# Load in FORMATTED TE PAV matrix for inbreds and hybrids
# te_genotypes_inbreds <-  data.table::fread("/workdir/mbb262/genotypes/hybrid_numeric_geno_no_header.txt", nThread = numThreads)
te_genotypes_inbreds <-  data.table::fread("/workdir/mbb262/genotypes/NAM_anchorwave_indels.PAV01.noAncestralSpecification.2022-09-27_rformat.txt", nThread = numThreads)
colnames(te_genotypes_inbreds)[1] <- "Taxa"

# Just subset to inbred parents for test
# te_genotypes_inbreds <- temp_te[1:26,]


# Load in expression -----------------------------------------------------------

# Load in expression data (taxa x genes)
expression_df <- data.table::fread("/workdir/mbb262/nam_hybrid_rnaseq/output/counts/nam_counts_vst.csv", nThread = numThreads)

# Format to hapmap IDs
x_ref_nam_hapmap <- read.csv("/home/mbb262/git_projects/te_ase_nam/data/nam_hapmap_xref.csv")
sampleTable <- read.csv("/workdir/mbb262/nam_hybrid_rnaseq/sample_metadata.csv") %>% 
  select(sample_name, cultivar, age, tissue)

# Change some sample names to mathc between files
x_ref_nam_hapmap$NAM <- gsub("B47", "PHB47", x_ref_nam_hapmap$NAM)
x_ref_nam_hapmap$NAM <- gsub("HP301", "Hp301", x_ref_nam_hapmap$NAM)
sampleTable <- merge(x_ref_nam_hapmap, sampleTable, by.x = "NAM", by.y = "cultivar") %>% 
  select(-NAM)

# Check dimensions before merge
dim(sampleTable)
dim(expression_df)

# Separate expression by tissue
expression_df <- base::merge(sampleTable, expression_df, by.y = "Taxa", by.x = "sample_name") 
dim(expression_df)

expression_gp <- expression_df %>% filter(tissue == "Growing point")
expression_lt <- expression_df %>% filter(tissue == "V3 leaf tip")
expression_lb <- expression_df %>% filter(tissue == "V3 leaf base")
expression_flf50 <- expression_df %>% filter(tissue == "Middle section of adult flag leaf" & age == "50 days")
expression_flf58 <- expression_df %>% filter(tissue == "Middle section of adult flag leaf" & age == "58 days")

# Check dimensions
dim(expression_gp)
dim(expression_lt)
dim(expression_lb)
dim(expression_flf50)
dim(expression_flf58)

# Peep the first part of the data
expression_gp[,1:5]
expression_lt[,1:5]
expression_lb[,1:5]
expression_flf50[,1:5]
expression_flf58[,1:5]

# Remove extra columns
expression_gp <- expression_gp %>% select(-sample_name, -age, -tissue, -sample_name)
expression_lt <- expression_lt %>% select(-sample_name, -age, -tissue, -sample_name)
expression_lb <- expression_lb %>% select(-sample_name, -age, -tissue, -sample_name)
expression_flf50 <- expression_flf50 %>% select(-sample_name, -age, -tissue, -sample_name)
expression_flf58 <- expression_flf58 %>% select(-sample_name, -age, -tissue, -sample_name)

# Change column names
colnames(expression_gp)[1] <- "Taxa"
colnames(expression_lt)[1] <- "Taxa"
colnames(expression_lb)[1] <- "Taxa"
colnames(expression_flf50)[1] <- "Taxa"
colnames(expression_flf58)[1] <- "Taxa"


# Create windows around genes for maize K matrices -----------------------------

# Download files
# system("wget -P /workdir/mbb262/ https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz")
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
te_key_file <- read.csv("/workdir/mbb262/genotypes/te_key_file.csv")
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

# Remove blacklisted hiccup genes
test_genes <- test_genes[!grepl('Zm00001eb027000|Zm00001eb126900|Zm00001eb409470|Zm00001eb229190|Zm00001eb241370', test_genes)]

# Debug section
# grep("Zm00001eb379090", test_genes)
# gene <-1254

# Subset to specific genes to manually parallellize this
# test_genes <- test_genes[1:4222]
test_genes <- test_genes[11868:length(test_genes)] # help out model 3 and analyze the last 800 genes

# Create place to store results
te_model_effects_results <- list()
te_model_stats_results <- list()

# Run loop; length(test_genes)
for (gene in 1:length(test_genes)) {
  
  # Subset out single gene's coordinates
  maize_coord <- gff_5utr %>% filter(attributes == test_genes[gene])
  print(test_genes[gene])
  
  # Subset vcf file
  gr_maize_k <- GenomicRanges::makeGRangesFromDataFrame(maize_coord)
  maize_k_vcf <- tasGeno %>% rTASSEL::filterGenotypeTableSites(siteRangeFilterType = "none", gRangesObj = gr_maize_k)
  
  # Calculate maize K
  k_maize <- kinshipMatrix(tasObj = maize_k_vcf)

  # Filter TE genotypes to those near the [gene] of interest
  gene_tes <- gene_te_overlaps %>% data.frame() %>% filter(v5_gene == test_genes[gene]) %>% select(te_marker_id)
  
  # If there are no TEs, move on to next gene, else do the models
  if (nrow(gene_tes) != 0){
    gene_tes <- c("Taxa", gene_tes$te_marker_id)
    filter_te_geno <- te_genotypes_inbreds[, ..gene_tes]
    dim(filter_te_geno)
    colnames(filter_te_geno)[1] <- "<Marker>"
    
    # filter_te_geno <- tibble::rownames_to_column(filter_te_geno, "<Marker>")
    data.table::fwrite(filter_te_geno, "/workdir/mbb262/01_temp_genotype.txt", sep = "\t")
    system("echo '<Numeric>' | cat - /workdir/mbb262/01_temp_genotype.txt > temp && mv temp /workdir/mbb262/01_temp_genotype.txt")
    
    # For each of the four tissues ---------------------------
    
    # Growing point ---------------
    # Filter expression matrix to just that gene
    filtered_expression_df_gp <- expression_gp %>% select(Taxa, test_genes[gene])
    
    # Turn genotypes and phenotypes into a tassel object
    tasGenoPheno_gp <- rTASSEL::readGenotypePhenotype(genoPathOrObj = "/workdir/mbb262/01_temp_genotype.txt",
                                                      phenoPathDFOrObj = filtered_expression_df_gp, 
                                                      taxaID = "Taxa", attributeTypes = NULL)
    
    # Run mixed model
    tasMLM_gp <- rTASSEL::assocModelFitter(tasObj = tasGenoPheno_gp, 
                                           formula = . ~ ., fitMarkers = TRUE, kinship = k_maize,
                                           fastAssociation = FALSE,maxThreads = numThreads)
    
    # Leaf tip --------------------
    filtered_expression_df_lt <- expression_lt %>% select(Taxa, test_genes[gene])
    tasGenoPheno_lt <- rTASSEL::readGenotypePhenotype(genoPathOrObj = "/workdir/mbb262/01_temp_genotype.txt",
                                                      phenoPathDFOrObj = filtered_expression_df_lt, 
                                                      taxaID = "Taxa", attributeTypes = NULL)
    tasMLM_lt <- rTASSEL::assocModelFitter(tasObj = tasGenoPheno_lt, 
                                           formula = . ~ ., fitMarkers = TRUE, kinship = k_maize,
                                           fastAssociation = FALSE,maxThreads = numThreads)
    
    # Leaf base -------------------
    filtered_expression_df_lb <- expression_lb %>% select(Taxa, test_genes[gene])
    tasGenoPheno_lb <- rTASSEL::readGenotypePhenotype(genoPathOrObj = "/workdir/mbb262/01_temp_genotype.txt",
                                                      phenoPathDFOrObj = filtered_expression_df_lb, 
                                                      taxaID = "Taxa", attributeTypes = NULL)
    tasMLM_lb <- rTASSEL::assocModelFitter(tasObj = tasGenoPheno_lb,
                                           formula = . ~ ., fitMarkers = TRUE, kinship = k_maize,
                                           fastAssociation = FALSE,maxThreads = numThreads)
    # Flag leaf at flowering ------
    filtered_expression_df_flf50 <- expression_flf50 %>% select(Taxa, test_genes[gene])
    tasGenoPheno_flf50 <- rTASSEL::readGenotypePhenotype(genoPathOrObj = "/workdir/mbb262/01_temp_genotype.txt",
                                                         phenoPathDFOrObj = filtered_expression_df_flf50, 
                                                         taxaID = "Taxa", attributeTypes = NULL)
    tasMLM_flf50 <- rTASSEL::assocModelFitter(tasObj = tasGenoPheno_flf50, 
                                              formula = . ~ ., fitMarkers = TRUE, kinship = k_maize,
                                              fastAssociation = FALSE,maxThreads = numThreads)
    
    # Add a tissue label to mixed model results --------------
    tasMLM_gp$MLM_Effects$tissue <- rep("Growing Point", nrow(tasMLM_gp$MLM_Effects))
    tasMLM_gp$MLM_Stats$tissue <- rep("Growing Point", nrow(tasMLM_gp$MLM_Stats))
    tasMLM_lt$MLM_Effects$tissue <- rep("Leaf Tip", nrow(tasMLM_lt$MLM_Effects))
    tasMLM_lt$MLM_Stats$tissue <- rep("Leaf Tip", nrow(tasMLM_lt$MLM_Stats))
    tasMLM_lb$MLM_Effects$tissue <- rep("Leaf Base", nrow(tasMLM_lb$MLM_Effects))
    tasMLM_lb$MLM_Stats$tissue <- rep("Leaf Base", nrow(tasMLM_lb$MLM_Stats))
    tasMLM_flf50$MLM_Effects$tissue <- rep("Flaf Leaf at Flowering", nrow(tasMLM_flf50$MLM_Effects))
    tasMLM_flf50$MLM_Stats$tissue <- rep("Flaf Leaf at Flowering", nrow(tasMLM_flf50$MLM_Stats))
    
    # Combine results across models
    te_model_effects_results[[gene]] <- data.table::as.data.table(rbind(tasMLM_gp$MLM_Effects, tasMLM_lt$MLM_Effects, tasMLM_lb$MLM_Effects, tasMLM_flf50$MLM_Effects))
    te_model_stats_results[[gene]] <- data.table::as.data.table(rbind(tasMLM_gp$MLM_Stats, tasMLM_lt$MLM_Stats, tasMLM_lb$MLM_Stats, tasMLM_flf50$MLM_Stats))
  }
  
  # insurance policy
  # if (gene %in% c(2, 500,1000, 1500, 2000, 2500, 3000, 3500, 4000)){
  #   data.table::fwrite(rbindlist(te_model_effects_results), paste0("/workdir/mbb262/nam_hybrid_rnaseq/output/model_results/01_insurance_save_mlm_effects_", gene, ".txt"), nThread = numThreads)
  #   data.table::fwrite(rbindlist(te_model_stats_results), paste0("/workdir/mbb262/nam_hybrid_rnaseq/output/model_results/01_insurance_save_mlm_stats_", gene, ".txt"), nThread = numThreads)
  # }
  if (gene %in% c(100, 500,750)){
    data.table::fwrite(rbindlist(te_model_effects_results), paste0("/workdir/mbb262/nam_hybrid_rnaseq/output/model_results/03_end_insurance_save_mlm_effects_", gene, ".txt"), nThread = numThreads)
    data.table::fwrite(rbindlist(te_model_stats_results), paste0("/workdir/mbb262/nam_hybrid_rnaseq/output/model_results/03_end_insurance_save_mlm_stats_", gene, ".txt"), nThread = numThreads)
  }
}

# Combine lists
combined_te_model_effects_results <- rbindlist(te_model_effects_results)
combined_te_model_stats_results <- rbindlist(te_model_stats_results)

# Write to file
# data.table::fwrite(combined_te_model_effects_results, "/workdir/mbb262/nam_hybrid_rnaseq/output/model_results/01_te_model_effects_results.csv", nThread = numThreads)
# 
# data.table::fwrite(combined_te_model_stats_results, "/workdir/mbb262/nam_hybrid_rnaseq/output/model_results/01_te_model_stats_results.csv", nThread = numThreads)
data.table::fwrite(combined_te_model_effects_results, "/workdir/mbb262/nam_hybrid_rnaseq/output/model_results/03_end_te_model_effects_results.csv", nThread = numThreads)

data.table::fwrite(combined_te_model_stats_results, "/workdir/mbb262/nam_hybrid_rnaseq/output/model_results/03_end_te_model_stats_results.csv", nThread = numThreads)


