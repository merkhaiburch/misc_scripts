# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch + Aimee Schutz
# Contact... mbb262@cornell.edu
# Date...... 2022-08-05 
# Updated... 2022-08-05

# Description 
# Prepare matrices for tassel, run tassel gui, plot results in R
# ---------------------------------------------------------------

# Load in kinship matrix, format, export
k <- read.csv("~/Downloads/germano/Original99_kmatrix.csv")
dim(k)
colnames(k) <- gsub("\\.", "_", colnames(k))
rownames(k) <- gsub("-", "_", rownames(k))
write.table(k, "~/Downloads/germano/k.txt", quote = F, row.names = T, sep = "\t", col.names = F)

# load in phenotypes, format, export
p <- read.csv("~/Downloads/germano/Original99_phenomatrix.csv")
p <- p[,c(1,10)]
colnames(p) <- c("Taxa", "pheno")
p$Taxa <- gsub("-", "_", p$Taxa)
write.table(p, "~/Downloads/germano/p.txt", quote = F, row.names = F, sep = "\t")

# Load in genotypes, format, export
x <- data.table::fread("~/Downloads/germano/Original99_pavMatrix.csv") %>% data.frame()
x$taxa <- gsub("-", "_", x$taxa)
x <- x[, grep("_T001|taxa", colnames(x))]

# Calculate row sums
temp <- data.frame(colnames(x[,-1]), colSums(x[,-1]))
hist(temp[,2])
colnames(temp) <- c("taxa", "freq")

# filter on a maf like metric
temp <- temp %>% filter(freq >= 6)
hist(temp$freq)
temp2 <- temp$taxa
temp2 <- c("taxa", temp2)

# Subset PAV matrix to sites with more than 6 taxa having that allele
lala <- x %>% select(temp2)

# Write to file
data.table::fwrite(lala, "~/Downloads/germano/x.txt", quote = F, row.names = F, sep = "\t")

# -----------------------------------------------------------------------
# Plot results
# -----------------------------------------------------------------------

# Load file into r from tassel gui
gwas_result <- read.delim("~/Desktop/TutorialData/perennial_mlm.txt")
gwas_result <- gwas_result[-1,]

# Load in gff file, parse out transcript names
gff <- ape::read.gff("~/Downloads/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3", 
                     na.strings = c(".", "?"), 
                     GFF3 = TRUE)
gff <- gff %>% filter(type == "mRNA")
gff$attributes <- gsub("ID=", "", gff$attributes)
gff$attributes <- gsub(";.*", "", gff$attributes)
gff <- gff %>% filter(seqid %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10"))

# Join gwas results and gff file for position information
gwas_result <- merge(gwas_result, gff, by.x = "Marker", by.y = "attributes", all.x = TRUE)

# Format x axis to plot by chromosome
data_cum <- gwas_result %>% 
  group_by(seqid) %>% 
  summarise(max_bp = max(start)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  select(seqid, bp_add)
gwas_result <- gwas_result %>% 
  inner_join(data_cum, by = "seqid") %>% 
  mutate(bp_cum = start + bp_add)

# filter crazy p-value outliers
gwas_result_filter <- gwas_result %>% filter(p > 1.8729e-115)

# Plot out
a <- ggplot(gwas_result_filter, aes(x = bp_cum, y = -log10(p), color = as.factor(seqid))) +
  geom_point(size = 1) +
  labs(x = NULL, y = "-log10(p)", title = "perennial ~ PAV (>6 taxa) + K") + 
  theme_bw()+
  theme(legend.position="bottom")
ggsave(plot = a,filename = "~/Downloads/germano/gwas_plot_aimee.png")
