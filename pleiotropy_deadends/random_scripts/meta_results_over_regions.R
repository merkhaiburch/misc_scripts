# ---------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2020-06-30 
#
# Description 
#   - Assign Plant Trait Ontologies onto phenotypes collected in:
#   - NAM, Goodman Panel, more TBA
# ---------------------------------------------------------------

# Load packages
library(dplyr)
library(ggplot2)
library(data.table)
library(RColorBrewer)

# # Load data
# gene <- read.csv("~/Downloads/gene_hits.csv", header = TRUE) %>% data.table::as.data.table()
# intergene <- read.csv("~/Downloads/intergene_hits.csv", header = TRUE)
# 
# # Sum traits by chromosome
# intergene <- intergene %>% 
#   group_by(trait) %>% 
#   summarise(snp_hits = sum(snp_hits)) %>% 
#   data.table::as.data.table()
# 
# # Change colnames before merge
# colnames(intergene)[2] <- "intergenic_hits"
# colnames(gene)[2] <- "genic_hits"
# 
# # Join by trait
# data.table::setkey(gene, trait)
# data.table::setkey(intergene, trait)
# joined_gi <- intergene[gene, nomatch = 0]
# 
# # Melt 
# melt_joined_gi <- reshape::melt(joined_gi, id = "trait")
# 
# # Plot
# ggplot(melt_joined_gi, aes(x = variable, y = log(value), fill = variable)) +
#   geom_boxplot() +
#   theme(legend.position="none") +
#   xlab("Genomic Region") +
#   ylab("log(Number of hits >1e-6)")
# 
# # Annotate factors
# annotate_traits <- melt_joined_gi %>% filter(!grepl("Zhou_2019", trait))
# annotate_metabolites <- melt_joined_gi %>% filter(grepl("Zhou_2019", trait))
# annotate_traits$trait_type <- paste0("classical")
# annotate_metabolites$trait_type <- paste0("metabolite")
# melt_joined_gi$trait_type <- paste0("all traits")
# 
# # Rbind all together
# all_melted_annotated <- rbind(melt_joined_gi, annotate_traits, annotate_metabolites)
# 
# # Plot
# ggplot(all_melted_annotated, aes(x = variable, y = log10(value), fill = variable)) +
#   geom_boxplot() +
#   theme(legend.position="none") +
#   xlab("Genomic Region") +
#   ylab("log10(Number of hits >1e-6)")
# 
# # Plot
# ggplot(all_melted_annotated, aes(x = trait_type, y = value, fill = variable))+
#   geom_bar(position="fill", stat="identity")
# 
# 
# # -------------------
# # Query gene regions
# # -------------------
# 
# # LOAD RESULTS
# utr_5 <- read.csv("~/Downloads/mgwas_results/five_utr_hits.csv", header = TRUE)
# utr_3 <- read.csv("~/Downloads/mgwas_results/three_utr_hits.csv", header = TRUE)
# exon <- read.csv("~/Downloads/mgwas_results/exon_hits.csv", header = TRUE)
# intron <- read.csv("~/Downloads/mgwas_results/intron_hits.csv", header = TRUE)
# 
# # Sum traits by chromosome
# utr_5 <- utr_5 %>% group_by(trait) %>% 
#   summarise(snp_hits = sum(snp_hits)) %>% data.table::as.data.table()
# utr_3 <- utr_3 %>% group_by(trait) %>% 
#   summarise(snp_hits = sum(snp_hits)) %>% data.table::as.data.table()
# exon <- exon %>% group_by(trait) %>% 
#   summarise(snp_hits = sum(snp_hits)) %>% data.table::as.data.table()
# intron <- intron %>% group_by(trait) %>% 
#   summarise(snp_hits = sum(snp_hits)) %>% data.table::as.data.table()
# 
# # Change colnames before merge
# colnames(utr_5)[2] <- "utr_5_hits"
# colnames(utr_3)[2] <- "utr_3_hits"
# colnames(exon)[2] <- "exon_hits"
# colnames(intron)[2] <- "intron_hits"
# 
# # Join by trait
# data.table::setkey(utr_5, trait)
# data.table::setkey(utr_3, trait)
# data.table::setkey(exon, trait)
# data.table::setkey(intron, trait)
# joined_genic_1 <- merge.data.table(utr_5, utr_3, all = TRUE)
# joined_genic_2 <- merge.data.table(exon, intron, all = TRUE)
# joined_genic <- merge.data.table(joined_genic_1, joined_genic_2, all = TRUE)
# 
# # Melt 
# melt_joined_genic <- reshape::melt(joined_genic, id = "trait")
# 
# # Add annotation for type of trait
# annotate_traits <- melt_joined_genic %>% filter(!grepl("Zhou_2019", trait))
# annotate_metabolites <- melt_joined_genic %>% filter(grepl("Zhou_2019", trait))
# annotate_traits$trait_type <- paste0("classical")
# annotate_metabolites$trait_type <- paste0("metabolite")
# melt_joined_genic$trait_type <- paste0("all traits")
# 
# # Rbind all together
# all_melted_annotated_genic <- rbind(melt_joined_genic, annotate_traits, annotate_metabolites) %>% na.omit()
# 
# # Plot
# ggplot(all_melted_annotated_genic, aes(x = trait_type, y = value, fill = variable))+
#   geom_bar(position="fill", stat="identity")
# 
# # Format NS data
# ns_exons <- read.csv("~/Downloads/mgwas_results/exon_hits_ns.csv") %>% 
#   summarise(snp_hits = sum(snp_hits)) %>% 
#   tibble::add_column(trait_type = paste0("input SNPs"), variable= paste0("exon_hits"))
# 
# ns_introns <- read.csv("~/Downloads/mgwas_results/intron_hits_ns.csv") %>% 
#   summarise(snp_hits = sum(snp_hits)) %>% 
#   tibble::add_column(trait_type = paste0("input SNPs"), variable= paste0("intron_hits"))
# 
# ns_intergenic <- read.csv("~/Downloads/mgwas_results/intergene_hits_ns.csv") %>% 
#   summarise(snp_hits = sum(snp_hits)) %>% 
#   tibble::add_column(trait_type = paste0("input SNPs"), variable= paste0("intergenic_hits"))
# 
# ns_utr_3 <- read.csv("~/Downloads/mgwas_results/three_utr_hits_ns.csv") %>% 
#   summarise(snp_hits = sum(snp_hits)) %>% 
#   tibble::add_column(trait_type = paste0("input SNPs"), variable= paste0("utr_3_hits"))
# 
# ns_utr_5 <- read.csv("~/Downloads/mgwas_results/five_utr_hits_ns.csv") %>% 
#   summarise(snp_hits = sum(snp_hits)) %>% 
#   tibble::add_column(trait_type = paste0("input SNPs"), variable= paste0("utr_5_hits"))
# 
# ns_results <- rbind(ns_exons, ns_introns, ns_intergenic, ns_utr_3, ns_utr_5)
# colnames(ns_results) <- c("value", "trait_type","variable")
# 
# # ----------------------
# # # All categories plot
# # ----------------------
# 
# # Combine all datasets with sig. SNPs with NS SNPs
# temp <- all_melted_annotated %>% filter(variable != "genic_hits") # replacing genic hits with introns/exons/utrs
# across_region_types <- rbind(temp, all_melted_annotated_genic, ns_results, fill = TRUE)
# 
# # Plot
# ggplot(across_region_types, aes(x = trait_type, y = value, fill = variable))+
#   geom_bar(position="fill", stat="identity")
# 
# # Get proportion by category
# temp <- across_region_types %>% 
#   group_by(variable, trait_type) %>% 
#   summarise(value = sum(value))
# 
# lala <- temp %>% 
#   group_by(trait_type) %>% 
#   mutate(freq = value / sum(value))
# 
# # Reorder factor levels
# lala$trait_type <- factor(lala$trait_type, levels = c("input SNPs", "all traits", "classical", "metabolite"))
# ggplot(lala, aes(x = trait_type, y = freq, fill = variable)) +
#   geom_bar(position = position_stack(), stat = "identity", width = .7) +
#   geom_text(aes(label = paste0(signif((freq*100),4), "%") ), position = position_stack(vjust = 0.5), size = 4) +
#   scale_fill_brewer(palette="Set3")
# 
# # Count total number of SNPs
# lala %>% group_by(trait_type) %>% summarise(value = sum(value))
# 
# 
# # -------------------------
# # Open vs closed chromatin
# # -------------------------
# 
# # Load files
# open_sig <- read.csv("~/Downloads/mgwas_results/open_chrom_hits.csv", header = TRUE)
# open_ns <- read.csv("~/Downloads/mgwas_results/open_chrom_hits_ns.csv", header = TRUE)
# closed_sig <- read.csv("~/Downloads/mgwas_results/closed_chrom_hits.csv", header = TRUE)
# closed_ns <- read.csv("~/Downloads/mgwas_results/closed_chrom_hits_ns.csv", header = TRUE)
# 
# # Sum results by trait across all chromosomes
# open_sig <- open_sig %>% group_by(trait) %>% 
#   summarise(snp_hits = sum(snp_hits)) %>% data.table::as.data.table()
# closed_sig <- closed_sig %>% group_by(trait) %>% 
#   summarise(snp_hits = sum(snp_hits)) %>% data.table::as.data.table()
# 
# # Rename columns
# colnames(open_sig)[2] <- "open_chrom_hits"
# colnames(closed_sig)[2] <- "closed_chrom_hits"
# 
# # Join
# data.table::setkey(open_sig, trait)
# data.table::setkey(closed_sig, trait)
# joined_chrom_status <- merge.data.table(open_sig, closed_sig, all = TRUE)
# 
# # Melt 
# melt_joined_chrom_status <- reshape::melt(joined_chrom_status, id = "trait")
# 
# # Add annotation for type of trait
# annotate_traits_chrom_status <- melt_joined_chrom_status %>% filter(!grepl("Zhou_2019", trait))
# annotate_metabolites_chrom_status <- melt_joined_chrom_status %>% filter(grepl("Zhou_2019", trait))
# annotate_traits_chrom_status$trait_type <- paste0("classical")
# annotate_metabolites_chrom_status$trait_type <- paste0("metabolite")
# melt_joined_chrom_status$trait_type <- paste0("all traits")
# 
# # Rbind all together
# all_melted_annotated_chrom_status <- rbind(melt_joined_chrom_status, annotate_traits_chrom_status, annotate_metabolites_chrom_status) %>% na.omit()
# 
# # Add in ns results (all input SNPs)
# closed_ns <- closed_ns %>% 
#   summarise(snp_hits = sum(snp_hits)) %>% 
#   tibble::add_column(trait_type = paste0("input SNPs"), variable= paste0("closed_chrom_hits")) %>% 
#   data.table::as.data.table()
# 
# open_ns <- open_ns %>% 
#   summarise(snp_hits = sum(snp_hits)) %>%
#   tibble::add_column(trait_type = paste0("input SNPs"), variable= paste0("open_chrom_hits")) %>% 
#   data.table::as.data.table()
# 
# colnames(closed_ns)[1] <- "value"
# colnames(open_ns)[1] <- "value"
# 
# # Get proportion by category
# temp <- all_melted_annotated_chrom_status %>% 
#   group_by(variable, trait_type) %>% 
#   summarise(value = sum(value)) %>% 
#   data.table::as.data.table()
# lala <- rbind(temp, open_ns, closed_ns) # Add in ns results
# lala <- lala %>% 
#   group_by(trait_type) %>% 
#   mutate(freq = value / sum(value))
# 
# # Reorder factor levels
# lala$trait_type <- factor(lala$trait_type, levels = c("input SNPs", "all traits", "classical", "metabolite"))
# ggplot(lala, aes(x = trait_type, y = freq, fill = variable)) +
#   geom_bar(position = position_stack(), stat = "identity", width = .7) +
#   geom_text(aes(label = paste0(signif((freq*100),4), "%") ), position = position_stack(vjust = 0.5), size = 4) +
#   scale_fill_brewer(palette="Set3")
# 
# 
# # ------------------------------------------
# # First attempt at pleiotropy distributions
# # ------------------------------------------
# 
# # Find files
# files_ply_ns <- list.files(path = "~/Documents/mgwas_results/", pattern = "pleio_hits_chrom*")
# files_ply_sig <- list.files(path = "~/Documents/mgwas_results/", pattern = "pleio_hits_stringent_chrom*")
# 
# # Read files in, give names
# ply_ns <- lapply(paste0("~/Documents/mgwas_results/", files_ply_ns), data.table::fread)
# ply_sig <- lapply(paste0("~/Documents/mgwas_results/", files_ply_sig), data.table::fread)
# names(ply_ns) <- files_ply_ns
# names(ply_sig) <- files_ply_sig
# 
# # Turn all into data.table objects
# ply_ns <- lapply(ply_ns, data.table::as.data.table)
# ply_sig <- lapply(ply_sig, data.table::as.data.table)
# 
# # Rbind all the lists together
# ply_ns <- data.table::rbindlist(ply_ns)
# ply_sig <- data.table::rbindlist(ply_sig)
# ply_ns$snp_type <- paste0("all_snps")
# ply_sig$snp_type <- paste0("sig_snps")
# ply <- rbind(ply_ns, ply_sig)
# 
# # Calculate median by chrom and type
# ply_stats <- ply %>% 
#   group_by(snp_seqid, snp_type) %>% 
#   summarise(num_trait = median(num_trait))
# 
# # Plot pleiotropy by chromosome facet
# ggplot(ply, aes(x=num_trait, color=snp_type)) +
#   geom_density(alpha=0.4, size = 1) +
#   facet_wrap(vars(snp_seqid), nrow = 5, ncol = 2) +
#   xlab("Number of traits associated with a single SNP") +
#   geom_vline(data=ply_stats, aes(xintercept=num_trait, color=snp_type),
#              linetype="dashed")
# 
# # Aggregate results across chromosomes
# ggplot(ply, aes(x=num_trait, color=snp_type)) +
#   geom_density(alpha=0.4, size = 1) +
#   xlab("Number of traits associated with a single SNP")
#   
# ply %>% 
#   group_by(snp_type) %>% 
#   summarise(num_trait = median(num_trait))
# 
# # return row with max value
# ply[which.max(ply$num_trait),]
# 
# # Plot only sig snps
# library(patchwork)
# temp1 <- ply %>% filter(!snp_type=="sig_snps")
# a <- ggplot(temp1, aes(x=num_trait)) +
#   geom_histogram(bins = 100) +
#   xlab("Number of traits associated with a single SNP")+
#   ylab("Number of SNPs") +
#   theme(legend.position = "none")
# 
# temp2 <- ply %>% filter(!snp_type=="sig_snps", num_trait > 0)
# b <- ggplot(temp2, aes(x=num_trait)) +
#   geom_histogram(bins = 100) +
#   xlab("Number of traits associated with a single SNP >0")+
#   ylab("Number of SNPs") +
#   theme(legend.position = "none")
# 
# a|b
# 
# 
# 
# ----------------------------------------------
# Results after filtering to top 1% of SNPs
# by haplotype, similar to peak calling
# ---------------------------------------------

# Load data
gene <- read.csv("~/Documents/mgwas_results/askdb_triplet_results/trip_gene_hits_hap_peaks.csv", header = TRUE) %>% data.table::as.data.table()
intergene <- read.csv("~/Documents/mgwas_results/askdb_triplet_results/trip_intergene_hits_hap_peaks.csv", header = TRUE)

# Sum traits by chromosome
intergene <- intergene %>%
  group_by(trait) %>%
  summarise(snp_hits = sum(snp_hits)) %>%
  data.table::as.data.table()
gene <- gene %>%
  group_by(trait) %>%
  summarise(snp_hits = sum(snp_hits)) %>%
  data.table::as.data.table()

# Change colnames before merge
colnames(intergene)[2] <- "intergenic_hits"
colnames(gene)[2] <- "genic_hits"

# Join by trait
data.table::setkey(gene, trait)
data.table::setkey(intergene, trait)
joined_gi <- intergene[gene, nomatch = 0]

# Melt
melt_joined_gi <- reshape::melt(joined_gi, id = "trait")

# Annotate factors
annotate_traits <- melt_joined_gi %>% filter(!grepl("Zhou_2019", trait))
annotate_metabolites <- melt_joined_gi %>% filter(grepl("Zhou_2019", trait))
annotate_traits$trait_type <- paste0("classical")
annotate_metabolites$trait_type <- paste0("metabolite")
melt_joined_gi$trait_type <- paste0("all traits")

# Rbind all together
all_melted_annotated <- rbind(melt_joined_gi, annotate_traits, annotate_metabolites)

# Plot
# ggplot(all_melted_annotated, aes(x = trait_type, y = value, fill = variable))+
#   geom_bar(position="fill", stat="identity")
#

# -------------------
# Query gene regions
# -------------------

# LOAD RESULTS
utr_5 <- read.csv("~/Documents/mgwas_results/askdb_triplet_results/trip_five_utr_hits_hap_peaks.csv", header = TRUE)
utr_3 <- read.csv("~/Documents/mgwas_results/askdb_triplet_results/trip_three_utr_hits_hap_peaks.csv", header = TRUE)
exon <- read.csv("~/Documents/mgwas_results/askdb_triplet_results/trip_exon_hits_hap_peaks.csv", header = TRUE)
intron <- read.csv("~/Documents/mgwas_results/askdb_triplet_results/trip_intron_hits_hap_peaks.csv", header = TRUE)

# Sum traits by chromosome
utr_5 <- utr_5 %>% group_by(trait) %>% 
  summarise(snp_hits = sum(snp_hits)) %>% data.table::as.data.table()
utr_3 <- utr_3 %>% group_by(trait) %>% 
  summarise(snp_hits = sum(snp_hits)) %>% data.table::as.data.table()
exon <- exon %>% group_by(trait) %>% 
  summarise(snp_hits = sum(snp_hits)) %>% data.table::as.data.table()
intron <- intron %>% group_by(trait) %>% 
  summarise(snp_hits = sum(snp_hits)) %>% data.table::as.data.table()

# Change colnames before merge
colnames(utr_5)[2] <- "utr_5_hits"
colnames(utr_3)[2] <- "utr_3_hits"
colnames(exon)[2] <- "exon_hits"
colnames(intron)[2] <- "intron_hits"

# Join by trait
data.table::setkey(utr_5, trait)
data.table::setkey(utr_3, trait)
data.table::setkey(exon, trait)
data.table::setkey(intron, trait)
joined_genic_1 <- merge.data.table(utr_5, utr_3, all = TRUE)
joined_genic_2 <- merge.data.table(exon, intron, all = TRUE)
joined_genic <- merge.data.table(joined_genic_1, joined_genic_2, all = TRUE)

# Melt 
melt_joined_genic <- reshape::melt(joined_genic, id = "trait")

# Add annotation for type of trait
annotate_traits <- melt_joined_genic %>% filter(!grepl("Zhou_2019", trait))
annotate_metabolites <- melt_joined_genic %>% filter(grepl("Zhou_2019", trait))
annotate_traits$trait_type <- paste0("classical")
annotate_metabolites$trait_type <- paste0("metabolite")
melt_joined_genic$trait_type <- paste0("all traits")

# Rbind all together
all_melted_annotated_genic <- rbind(melt_joined_genic, annotate_traits, annotate_metabolites) %>% na.omit()

# Remove ionomics data
all_melted_annotated_genic <- all_melted_annotated_genic %>% filter(!grepl("IonomicsNamCombined", trait))

# Plot
ggplot(all_melted_annotated_genic, aes(x = trait_type, y = value, fill = variable))+
  geom_bar(position="fill", stat="identity")

# Format all SNPs in the database (labeled as ns  = non-sig but some are sig)
ns_exons <- read.csv("~/Documents/mgwas_results/exon_hits_ns.csv") %>% 
  summarise(snp_hits = sum(snp_hits)) %>% 
  tibble::add_column(trait_type = paste0("input SNPs"), variable= paste0("exon_hits"))

ns_introns <- read.csv("~/Documents/mgwas_results/intron_hits_ns.csv") %>% 
  summarise(snp_hits = sum(snp_hits)) %>% 
  tibble::add_column(trait_type = paste0("input SNPs"), variable= paste0("intron_hits"))

ns_intergenic <- read.csv("~/Documents/mgwas_results/intergene_hits_ns.csv") %>% 
  summarise(snp_hits = sum(snp_hits)) %>% 
  tibble::add_column(trait_type = paste0("input SNPs"), variable= paste0("intergenic_hits"))

ns_utr_3 <- read.csv("~/Documents/mgwas_results/three_utr_hits_ns.csv") %>% 
  summarise(snp_hits = sum(snp_hits)) %>% 
  tibble::add_column(trait_type = paste0("input SNPs"), variable= paste0("utr_3_hits"))

ns_utr_5 <- read.csv("~/Documents/mgwas_results/five_utr_hits_ns.csv") %>% 
  summarise(snp_hits = sum(snp_hits)) %>% 
  tibble::add_column(trait_type = paste0("input SNPs"), variable= paste0("utr_5_hits"))

ns_results <- rbind(ns_exons, ns_introns, ns_intergenic, ns_utr_3, ns_utr_5)
colnames(ns_results) <- c("value", "trait_type","variable")


# ----------------------
# # All categories plot
# ----------------------

# Combine all datasets with sig. SNPs with NS SNPs
temp <- all_melted_annotated %>% filter(variable != "genic_hits") # replacing genic hits with introns/exons/utrs
across_region_types <- rbind(temp, all_melted_annotated_genic, ns_results, fill = TRUE)

# Plot
ggplot(across_region_types, aes(x = trait_type, y = value, fill = variable))+
  geom_bar(position="fill", stat="identity")

# Get proportion by category
temp <- across_region_types %>% 
  group_by(variable, trait_type) %>% 
  summarise(value = sum(value))

lala <- temp %>% 
  group_by(trait_type) %>% 
  mutate(freq = value / sum(value))

# Reorder factor levels
lala$trait_type <- factor(lala$trait_type, levels = c("input SNPs", "all traits", "classical", "metabolite"))
por <- ggplot(lala, aes(x = trait_type, y = freq, fill = variable)) +
  geom_bar(position = position_stack(), stat = "identity", width = .7) +
  geom_text(aes(label = paste0(signif((freq*100),4), "%") ), position = position_stack(vjust = 0.5), size = 2) +
  scale_fill_brewer(palette="Set3") +
  theme(text = element_text(size = 7), legend.position = "top", legend.text=element_text(size=4.5)) +
  labs(y = "Proportion", x = "SNP or Trait Category") 

# Export plot
ggsave(por, filename = "git_projects/haplotype_v_snp_gwas/images/proportion_gwas_hits_2.png",
       width = 4, height = 3, units = "in")

# Count total number of SNPs
lala %>% group_by(trait_type) %>% summarise(value = sum(value))


# -------------------
# Pleiotropy distributions
# -------------------------

# Find files
files_ply_ns <- list.files(path = "~/Documents/mgwas_results/", pattern = "pleio_hits_chrom*")
files_ply_sig <- list.files(path = "~/Documents/mgwas_results/askdb_triplet_results/", pattern = "pleio_hits_chrom*")

# Read files in, give names
ply_ns <- lapply(paste0("~/Documents/mgwas_results/", files_ply_ns), data.table::fread)
ply_sig <- lapply(paste0("~/Documents/mgwas_results/askdb_triplet_results/", files_ply_sig), data.table::fread)
names(ply_ns) <- files_ply_ns
names(ply_sig) <- files_ply_sig

# Turn all into data.table objects
ply_ns <- lapply(ply_ns, data.table::as.data.table)
ply_sig <- lapply(ply_sig, data.table::as.data.table)

# Rbind all the lists together
ply_ns <- data.table::rbindlist(ply_ns)
ply_sig <- data.table::rbindlist(ply_sig)
ply_ns$snp_type <- paste0("all_snps")
ply_sig$snp_type <- paste0("sig_snps")
ply <- rbind(ply_ns, ply_sig)

# Calculate median by chrom and type
ply_stats <- ply %>% 
  group_by(snp_seqid, snp_type) %>% 
  summarise(num_trait = median(num_trait))

# genome wide snp pleiotropy
ply_ns %>% 
  group_by(snp_type) %>% 
  summarise(num_trait = median(num_trait))
ply_ns %>% filter(num_trait > 0) %>% group_by(snp_type) %>% 
  summarise(num_trait = meadian(num_trait))
dim(ply_ns %>% filter(num_trait > 0))

# return row with max value
ply[which.max(ply$num_trait),]

# Plot only non-sig snps
library(patchwork)
# temp <- ply_ns[1:1000000,]
a <- ggplot(ply_ns, aes(x=num_trait, color = snp_type, fill = snp_type)) +
  geom_histogram(bins = 100, position="identity", alpha=0.5, color = "darkolivegreen4", fill = "darkolivegreen4") +
  # geom_vline(aes(xintercept=median(num_trait)), color="black", linetype="solid", size=1) +
  xlab("#traits associated with a single SNP")+
  ylab("Number of SNPs") +
  theme(legend.position = "none", text = element_text(size = 25))

temp2 <- ply_ns %>% filter(num_trait > 0) 
b <- ggplot(temp2, aes(x=num_trait, color = snp_type, fill = snp_type)) +
  geom_histogram(bins = 100, position="identity", alpha=0.5, color = "darkolivegreen4", fill = "darkolivegreen4") +
  geom_vline(aes(xintercept=median(num_trait)), color="black", linetype="solid", size=1) +
  xlab("#traits associated with a single SNP >0")+
  ylab("Number of SNPs") +
  theme(legend.position = "none", text = element_text(size = 25))

# a|b

ggsave(a, filename = "git_projects/haplotype_v_snp_gwas/images/genome_wide_gwas_hits_1.png",
       width = 4, height = 3, units = "in")
# ggsave(b, filename = "git_projects/haplotype_v_snp_gwas/images/genome_wide_gwas_hits_2.png", 
#        width = 4, height = 3, units = "in")
