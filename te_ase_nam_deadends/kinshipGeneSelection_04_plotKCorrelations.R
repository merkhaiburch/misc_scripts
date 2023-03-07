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
options(java.parameters = "-Xmx70g")
options(java.parameters = "-Xmx70g")
options(java.parameters = "-Xmx70g")
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


# Load in data -----------------------------------------------------------------

# Load in expression data (taxa x genes)


# Load in PAV matrix (taxa x TE PAV)


# Load in genotype table for inbreds and hybrids
# of SNPs: 10,156,186; of taxa: 113
# mkdir /workdir/mbb262/genotypes
# scp -r mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/mbb262/te_ase_nam/genotypes_for_k/v5_inbred_hybrid_genotypes/hmp321_282_agpv5_inbred_hybrid_all_chroms.vcf.gz /workdir/mbb262/genotypes
tasGeno <-  rTASSEL::readGenotypeTableFromPath("/workdir/mbb262/genotypes/hmp321_282_agpv5_inbred_hybrid_all_chroms.vcf.gz", keepDepth = FALSE)
print(tasGeno)


# Create windows around genes for maize K matrices -----------------------------

# Download gff files
# mkdir /workdir/mbb262/k_files
# wget -P /workdir/mbb262/k_files/ https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz
# wget -P /workdir/mbb262/k_files/ https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.canonical_transcripts
# wget -P /workdir/mbb262/k_files/ https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.canonical_transcripts
# wget -P /workdir/mbb262/k_files/ https://de.cyverse.org/anon-files//iplant/home/shared/NAM/NAM_Canu1.8_TE_annotation_03032022/B73.PLATINUM.pseudomolecules-v1.fasta.mod.EDTA.TEanno.split.gff3.gz


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

# Save to disk
write.table(gff_5utr[,-4], "/workdir/mbb262/k_files/gene_coords_500bp.txt",
            quote = F, row.names = F, sep = "\t", col.names = F)


# Create windows around genes for TE K matrices --------------------------------

# Load in TE gff file and filter
gff_te <- ape::read.gff("/workdir/mbb262/k_files/B73.PLATINUM.pseudomolecules-v1.fasta.mod.EDTA.TEanno.split.gff3.gz", na.strings = c(".", "?"), GFF3 = TRUE) %>% 
  filter(seqid %in% c("B73_chr1","B73_chr2","B73_chr3","B73_chr4","B73_chr5","B73_chr6","B73_chr7","B73_chr8","B73_chr9","B73_chr10")) %>% 
  select("seqid", "start", "end", "strand", "attributes")
gff_te$seqid <- gsub("B73_", "", gff_te$seqid)

# Format gene gff file
gff_gene <- ape::read.gff("/workdir/mbb262/k_files/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz", 
                          na.strings = c(".", "?"), GFF3 = TRUE) %>% 
  filter(type == "gene", seqid %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10")) %>% 
  select("seqid", "start", "end", "strand", "attributes")
gff_gene$attributes <- gsub("ID=", "", gff_gene$attributes) # parse out gene names
gff_gene$attributes <- gsub(";.*", "", gff_gene$attributes) # parse out gene names

# Add a window to gene's start
te_window <- 5000
gff_gene$start <- gff_gene$start - te_window

# Find overlapping TEs within a "te_window" size of a gene
gff_gene <- data.table::as.data.table(gff_gene)
gff_te <- data.table::as.data.table(gff_te)
setkey(gff_gene, seqid, start, end)
gene_te_overlaps <- data.table::foverlaps(gff_te, gff_gene, type="within",nomatch=0L) #0L drops TEs with no overlapping genes

# Select just the TE columns
gene_te_overlaps <- gene_te_overlaps %>% select("seqid", "i.start", "i.end", "i.strand", "attributes", "i.attributes")

# Fix column names
colnames(gene_te_overlaps) <- c("seqid", "start", "end", "strand", "v5_gene", "attributes")

# Write TE file to disk
write.table(gene_te_overlaps[,-c(4,6)], "/workdir/mbb262/k_files/te_gene_overlaps.txt", 
            quote = F, row.names = F, sep = "\t", col.names = F)


# Correlate K matrices ---------------------------------------------------------

##################################
# Note: this is just a test, for the real, fast version, see kotlin script
##################################

# Example correlating lower triangle of matrices, removing the diagonal 
# (m2 <- matrix(1:5, 5, 5))
# (m2 <- m2[lower.tri(m2, diag = F)])
# (m1 <- matrix(-6.3:-10.3, 5, 5))
# (m1 <- m1[lower.tri(m1, diag = F)])
# cor(m1,m2)

# Get intersecting genes to iterate through
all_genes <- intersect(gene_te_overlaps$v5_gene, gff_5utr$attributes)

# Write to file
temp <- data.frame(all_genes)
# write.table(temp, "/workdir/mbb262/k_files/all_genes.txt", quote = F, row.names = F, col.names = F)

# Store results
k_correlations <- c()

# Run parallel loop
start_time_foreach <- Sys.time()
for (gene in 1:length(all_genes)) {
  # K maize functional -------------------
  # Subset out gene's coordinates
  maize_coord <- gff_5utr %>% filter(attributes == all_genes[gene])
  print(all_genes[gene])
  # Subset vcf file
  gr_maize_k <- GenomicRanges::makeGRangesFromDataFrame(maize_coord)
  maize_k_vcf <- tasGeno %>% rTASSEL::filterGenotypeTableSites(siteRangeFilterType = "none", gRangesObj = gr_maize_k)
  
  # if there are no sites, skip this gene
  if(!is.na(maize_k_vcf)){
    # Calculate maize K
    k_maize <- kinshipMatrix(tasObj = maize_k_vcf)
    k_maize_matrix <- as.matrix(k_maize)
    
    # Save number of SNPs
    maize_num_snps <- maize_k_vcf@jPositionList$numberOfSites()
    
    # Calculate K TE ----------------------
    # Subset out gene's coordinates
    te_coord <- gene_te_overlaps %>% filter(v5_gene == all_genes[gene])
    num_tes <- nrow(te_coord)
    
    # subset vcf file
    gr_te_k <- GenomicRanges::makeGRangesFromDataFrame(te_coord)
    te_k_vcf <- tasGeno %>% rTASSEL::filterGenotypeTableSites(siteRangeFilterType = "none", gRangesObj = gr_te_k)
    
    if (!is.na(te_k_vcf)){
      # Calculate TE K
      k_te <- kinshipMatrix(tasObj = te_k_vcf)
      k_te_matrix <- as.matrix(k_te)
      
      # Save number of SNPs and number of TEs
      te_num_snps <- te_k_vcf@jPositionList$numberOfSites()
      
      # Gather lower diagnoal of matrix
      kin1 <- k_maize_matrix[lower.tri(k_maize_matrix, diag = F)]
      kin2 <- k_te_matrix[lower.tri(k_te_matrix, diag = F)]
      
      # Correlate two K matrices
      cor_out <- cor(kin1, kin2)
      
      # Save output
      k_correlations <- rbind(k_correlations, data.frame(all_genes[gene], cor_out, maize_num_snps, te_num_snps, num_tes))
    } 
  }
}

# Estimate timing
end_time_foreach <- Sys.time()
(end_time_foreach - start_time_foreach)

# Write to file
write.csv(k_correlations, "/workdir/mbb262/k_correlations.csv", row.names = F, quote = F)


# Plotting ---------------------------------------------------------------------

# Read in file
k_correlations <- read.csv("/workdir/mbb262/k_correlations_merritt.csv")
k_correlations <- read.delim("~/Downloads/kinship_correlation_500bpgene_5kbte.txt")
# k_correlations <- read.csv("/workdir/mbb262/kinship_correlation_terry.txt", sep = "\t")
# colnames(k_correlations)[2] <- "cor_out"

# Make R2
k_correlations$cor_out <- k_correlations$Pearson.Correlation^2

# Calculate mean
mu <- k_correlations %>% summarize(mean = mean(cor_out))

# Make plot
# title = expression(paste("Pearson ", R^{2}, " correlation between 500 bp maize regulatory and 5 Kb TE regulatory regions"))
library(ggplot2)
axis_text_size <- 22
plot1 <- ggplot(k_correlations, aes(x = cor_out)) +
  geom_histogram(position="identity", alpha=0.5, bins = 75)+
  geom_vline(data=mu, aes(xintercept=mean), linetype="dashed")+
  geom_vline(data=mu, aes(xintercept=0.9), linetype="solid")+
  labs(x = expression(paste("Pearson ", R^{2}, " correlation")),
       y = "Count") +
  theme_bw() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size)) 

# save
ggsave("/workdir/mbb262/pearson_cor_k.png", plot1, height = 10, width = 11.5, units = "in")

# count the number of genes before 0.9 threshold
k_correlations %>% select(cor_out) %>% filter(cor_out < 0.9) %>% summarise(n= n())


# See how correlations relate with other metrics -------------------------------

# ALl by all correlation matrix of values
PerformanceAnalytics::chart.Correlation(k_correlations[,-1], histogram=TRUE, pch=19)


# Compare with terry's version -------------------------------------------------

# Load in terry's file
terry <- read.delim("/workdir/mbb262/kinship_correlation_terry.txt")

# dim
dim(k_correlations)
dim(terry)

# head
head(k_correlations)
head(terry)

# merge
mer_terry <- merge(k_correlations, terry, by.x = "all_genes.gene.", by.y = "Gene.Name")
head(mer_terry)

# Correlate number of maize SNPs
cor(mer_terry$maize_num_snps, mer_terry$Number.of.SNPs.First)
plot(mer_terry$maize_num_snps, mer_terry$Number.of.SNPs.First)
table(mer_terry$maize_num_snps == mer_terry$Number.of.SNPs.First)

# Correlate number of TE SNPs
cor(mer_terry$te_num_snps, mer_terry$Number.of.SNPs.Second)
plot(mer_terry$te_num_snps, mer_terry$Number.of.SNPs.Second)
table(mer_terry$te_num_snps == mer_terry$Number.of.SNPs.Second)

# Correlate number of TEs
cor(mer_terry$num_tes, mer_terry$Number.Entries.in.Second.Gene.Location.File)
plot(mer_terry$num_tes, mer_terry$Number.Entries.in.Second.Gene.Location.File)
table(mer_terry$num_tes == mer_terry$Number.Entries.in.Second.Gene.Location.File)

# Correlate correlations
cor(mer_terry$cor_out, mer_terry$Pearson.Correlation)
plot(mer_terry$cor_out, mer_terry$Pearson.Correlation)
table(mer_terry$cor_out == mer_terry$Pearson.Correlation)

temp <- mer_terry %>% select(cor_out, Pearson.Correlation)

# They're all exactly the same and correlated at R = 1!
