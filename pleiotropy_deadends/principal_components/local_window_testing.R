# --------------------------------------------
# Merritt Burch
# mbb262@cornell.edu
# 2019-11-22
#
# Script to do local PCA in Ames with
# - fixed window sizes
# - variable number of genes/window
# --------------------------------------------

# Save image
# save.image("~/Box\ Sync/Cornell_PhD/labProjects/hap_gwa/ames_tests/2019_10_30_ames2NAM_loading_test/2019_11_15_605.Rdata")

# Load in source scripts
source('~/git_projects/pleiotropy/src/R/03_ames2any_matrix.R')
source('~/git_projects/pleiotropy/src/R/04_local_window_funs.R')

# Load packages
library(gdsfmt)
library(SNPRelate)
library(ggplot2)
library(ape)
library(dplyr)
library(broom)
library(magrittr)

# ------------------
#   Load test data
# ------------------

# Load in Beagle imputed SNPs that were filtered with the 'temp'
vcf.fn <- "~/Box Sync/Cornell_PhD/labProjects/hackathons/2019-10-21_hackathon/ZeaGBSv27_20171204_AGPv4_Ames_filteredGenotypes_beagleImputation_Filter_intersectSitesWithNAM_Filter.vcf"
beagle_vcf <- snpgdsVCF2GDS(vcf.fn, "beagle_filtered.gds", method="copy.num.of.ref", snpfirstdim=TRUE, ignore.chr.prefix = "S")
snpgdsSummary("beagle_filtered.gds")
genofile_ames <- snpgdsOpen(beagle_vcf)

# Load in NAM SNPs that are shared (intersect) with the Ames SNPs (did intersection in TASSEL)
nam.vcf <- "~/Box Sync/Cornell_PhD/labProjects/hackathons/2019-10-21_hackathon/intersectJoin_NAMbyAmes_shared_beagleImputed.vcf.gz"
nam_vcf <- snpgdsVCF2GDS(nam.vcf, "nam.gds", method="copy.num.of.ref", snpfirstdim=TRUE, ignore.chr.prefix = "S")
snpgdsSummary("nam.gds")
genofile_nam <- snpgdsOpen(nam_vcf)


# ------------------------------------------
#     Calculate ames2nam global PCs
#     uses ames2any_matrix functions
# ------------------------------------------

# Calculate Ames global PCs, get loadings
global_ames_BV <- get_ames_BV(X_ames = genofile_ames, Q = NULL, num_PCs = 3)

# Combine genofile objects
X_all <- combine_genofile_matrix(genofile_ames = genofile_ames, genofile_any = genofile_nam)

# Transfer loadings and coefficients over to all X (Ames, NAM) to get adjusted global PCs
ames2nam_gPCs <- ames2any(X_all = X_all$X_all,
                          Q = NULL,
                          B = global_ames_BV$B, 
                          V = global_ames_BV$V)

# Get Q matrix (population structure)
Q_all <- cbind(1, ames2nam_gPCs)


# ----------------------------------------
#      Calculate nam2nam global PCs
#         (no Ames adjustment)
# ----------------------------------------

# Calculate NAM global B and V matrix
global_nam_BV <- get_ames_BV(X_ames = genofile_nam, Q = NULL, num_PCs = 3)

# Make NAM X matrix and do formatting
X_all_nam_geno <- snpgdsGetGeno(genofile_nam, snpfirstdim = FALSE, with.id = TRUE, verbose = FALSE)
X_all_nam <- 2-X_all_nam_geno$genotype
colnames(X_all_nam) <- X_all_nam_geno$snp.id
rownames(X_all_nam) <- c(X_all_nam_geno$sample.id)
X_all_nam <- list(X_all = X_all_nam, snp_info = snpgdsSNPList(genofile_nam))

# Calculate PCs
nam2nam_gPCs <- ames2any(X_all = X_all_nam$X_all,
                         Q = NULL,
                         B = global_nam_BV$B, 
                         V = global_nam_BV$V)

# Make Q matrix (population structure)
Q_all_nam <- cbind(1, nam2nam_gPCs)


# ------------------------------------------
#    Calculate local PCs in fixed windows
# ------------------------------------------


# Use fixed widow function on ames to nam calculated PCs
fixed_ames2nam_local <- calc_local_pcs_fixed_window(X_all = X_all,
                                                    genofileID_ames = genofile_ames,
                                                    Q_all = Q_all,
                                                    numLocalPcs = 3,
                                                    windowSize = 1e7)

# Use fixed widow function on nam-only calculated PCs
fixed_nam2nam_local <- calc_local_pcs_fixed_window(X_all = X_all_nam,
                                                   genofileID_ames = genofile_nam,
                                                   Q_all = Q_all_nam,
                                                   numLocalPcs = 3,
                                                   windowSize = 1e7)


# ------------------------------------------
#   Calculate local PCs in defined windows
#           formatting step
# ------------------------------------------

# Get all maize v4 genes from gff file 
maize_gff <- read.gff("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/ames_tests/Sep25_localPCs/Zea_mays.B73_RefGen_v4.45.gff3" , na.strings = c(".", "?"), GFF3 = TRUE)

# Do some formatting to isolate gene names
maize_gff <- maize_gff[which(maize_gff$type == "gene"),]
maize_gff$gene <- gsub("ID=gene:", "", maize_gff$attributes)
maize_gff$gene <- gsub(";.*", "", maize_gff$gene)

# Rearrange, get rid of columns I don't need
genes <- maize_gff[, c(10,1,4,5)]
genes <- genes %>% 
  dplyr::filter(seqid %in% 1:10) %>% 
  dplyr::mutate(seqid = factor(.$seqid, levels = 1:10))

# Use the function to cut genome into 150 gene windows
# 177 genes = 226 windows; 178 = 225; 179=221; 180=220.33; 181 = 220!
window_file <- getWindowCoords(181, genes)

# Format file
window_file <- as.data.frame(window_file)
colnames(window_file) <- c("window", "chrom", "start", "stop")
window_file$chrom <- as.numeric(as.character(window_file$chrom))
window_file$start <- as.numeric(as.character(window_file$start))
window_file$stop <- as.numeric(as.character(window_file$stop))


# -------------------------------------------
# Local PC function with gene number windows
# -------------------------------------------

# Use on variable gene number windows with ames to nam PCs
geneWindow_ames2nam_local <- calc_local_pcs_gene_window(X_all = X_all,
                                                        genofileID_ames = genofile_ames,
                                                        Q_all = Q_all,
                                                        numLocalPcs = 3,
                                                        windowFile = window_file)

# Use on variable gene number windows with ames to nam PCs
geneWindow_nam2nam_local <- calc_local_pcs_gene_window(X_all = X_all_nam,
                                                       genofileID_ames = genofile_nam,
                                                       Q_all = Q_all_nam,
                                                       numLocalPcs = 3,
                                                       windowFile = window_file)


# ---------------------------------------------
# Calculate type III ss from these lPCs for NAM
# ---------------------------------------------

typeIIIss <- function(rowNames, globalPCs, localPCs){
  # Import file with all phenotypic data
  all_NAM_phenos <- read.csv("~/Box Sync/Cornell_PhD/labProjects/hackathons/2019-10-21_hackathon/all_NAM_phenos.txt", sep="")
  
  # Filter down phenotypes to only flowering time
  all_NAM_phenos <- all_NAM_phenos[,c(2,8:10)]

  # Combine global and local PCs, format taxa names
  nam_global_local <- cbind(rowNames, 
                            data.frame(globalPCs), 
                            data.frame(localPCs))
  
  # Change first column name
  colnames(nam_global_local)[1] <- "taxa"
  
  # Removes sequence identifiers from rownames
  nam_global_local[,1] <- gsub(":[0-9]{9}", "", nam_global_local[,1])
  
  # Merge dataframes together
  merged_pcs_phenos_nam <- merge(x = all_NAM_phenos, y = nam_global_local,
                                 by.x = "Geno_Code", by.y = "taxa")

  # Set options/contrasts that sets sums to zero
  options(contrasts = c("contr.sum", "contr.poly"))
  
  # Store the model:
  model_nam <- lm(Days_To_Silk_BLUP_Sum0607_Buckler2009 ~ ., merged_pcs_phenos_nam[,-c(1,2,4)])
  print(summary(model_nam))
  # The results give the type III SS, including the p-values from an F-test
  typeIII <- drop1(model_nam, .~., test = "F")
  
  # subset out 3 global PCs and intercept
  typeIII <- typeIII[-c(1:4),]
  
  # Make rownames one column
  typeIII <- data.table::setDT(typeIII, keep.rownames = TRUE)[]
  
  # Grep out chromosome, window, and PC info for plotting
  typeIII$chr <- gsub("_window_.*", "", typeIII$rn)
  typeIII$window <- gsub(".*window_", "", typeIII$rn)
  typeIII$PC <- gsub(".*_PC", "", typeIII$window)
  typeIII$window <- as.numeric(gsub("_PC[0-9]", "", typeIII$window))
  
  # Return the model, formatted results
  return(list(model_nam, typeIII))
}


# Fixed window methods
typeIII_fixed_ames2nam_local <- typeIIIss(rowNames = rownames(X_all$X_all),
                                          globalPCs = ames2nam_gPCs,
                                          localPCs = fixed_ames2nam_local)

typeIII_fixed_nam2nam_local <- typeIIIss(rowNames = rownames(X_all_nam$X_all),
                                         globalPCs = nam2nam_gPCs,
                                         localPCs = fixed_nam2nam_local)

# Gene window methods
typeIII_geneWindow_ames2nam_local <- typeIIIss(rowNames = rownames(X_all$X_all),
                                          globalPCs = ames2nam_gPCs,
                                          localPCs = geneWindow_ames2nam_local)

typeIII_geneWindow_nam2nam_local <- typeIIIss(rowNames = rownames(X_all_nam$X_all),
                                         globalPCs = nam2nam_gPCs,
                                         localPCs = geneWindow_nam2nam_local)


# Save all gene window PCs
tmp <- cbind(rownames(X_all_nam), ames2nam_gPCs, geneWindow_ames2nam_local)
write.table()



# ------------------------------------------
# Merge in flowering time loci for plotting
# ------------------------------------------

# Get data
dong_genes_with_coordinates <- read.csv("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/ames_tests/Sep25_localPCs/dong_genes_with_coordinates.csv", header = TRUE)

# More ugly code to find which fixed windows the dong genes are in for plotting
# Get window information
# Function gets coordiantes of window start and stops
get_windows <- function(genofileID, windowSize = 1e7){
  # Iterate across all chromosomes
  all_chr_formatted <- c()
  for (chrom in seq(1:10)){
    
    # Get whole SNP list
    snp_map <- snpgdsSNPList(genofileID)
    
    # Subset SNPs by chromosome and window
    windows_genome <- cut_width(snp_map$position[snp_map$chromosome == as.character(1)], width = 1e7)
    
    # Place in own dataframe
    temp <- data.frame(unique(windows_genome))
    temp <- stringr::str_split_fixed(temp$unique.windows_genome., ",", 2)
    
    # Add in chromosome
    temp <- data.frame(cbind(rep(chrom, nrow(temp)), seq(1:nrow(temp)), temp))
    
    # Remove brackets and other stuff
    temp[,3] <- gsub("\\(", "", temp[,3])
    temp[,3] <- gsub("\\[", "", temp[,3])
    temp[,4] <- gsub("]", "", temp[,4])
    
    # Into numeric
    temp[,3] <- as.numeric(as.character(temp[,3]))
    temp[,4] <- as.numeric(as.character(temp[,4]))
    
    # Add to collector
    all_chr_formatted <- rbind(all_chr_formatted, temp)
  }
  # Format collector
  all_chr_formatted <- data.frame(all_chr_formatted)
  return(all_chr_formatted)
}

# Use the function
ames_windows <- get_windows(genofile_ames, windowSize = 1e7)
colnames(ames_windows) <- c("chr", "window", "window_start", "window_stop")

# Subset SNPs by chromosome, window, put them into separate lists to iterate over
ranges <- merge(ames_windows, dong_genes_with_coordinates,by.x="chr", by.y = "seqid")
intersect_dong_ranges <- ranges[with(ranges, window_start <= start & window_stop >= end),]
intersect_dong_ranges$combined <- paste0("chr_", intersect_dong_ranges$chr, "_window_", intersect_dong_ranges$window, sep = "") # make column with info combined
intersect_dong_ranges$pos <- rep(150, nrow(intersect_dong_ranges))
intersect_dong_ranges <- intersect_dong_ranges[,-9] # Remove second chromosome column

# UGLY Code to get which of my variable number gene windows has Dong 2012 FT genes in (not elegant but works)
temp <- typeIII_geneWindow_ames2nam_local # copy data frame 
temp[[2]]$rn_noPC_ID <- gsub("_PC[0-9]", "", temp[[2]]$rn)
temp[[2]]$chromosome <- gsub("chr_", "", temp[[2]]$chr)
temp[[2]]$chromosome <- as.numeric(temp[[2]]$chromosome)
temp_dong_FTgenes <- merge(x = temp[[2]], y = window_file,by.x = "rn_noPC_ID", by.y = "window")
ranges <- merge(temp_dong_FTgenes, dong_genes_with_coordinates,by.x="chromosome", by.y = "seqid", allow.cartesian=TRUE) # Big and messy
intersect_dong_ranges <- ranges[with(ranges, start.x <= start.y & stop >= end),]
library(tidyverse)
unique_gene <- intersect_dong_ranges %>% distinct() # only keep distinct window/gene combinations
unique_gene <- unique_gene[,c(2,3,10:25)] # clean up
unique_gene$pos <- rep(150, nrow(unique_gene))

# Fixed window
typeIII_fixed_ames2nam_local[[2]]$rn2 <- gsub("_PC[0-9]", "", typeIII_fixed_ames2nam_local[[2]]$rn)
fixed_ames2nam_dong <- merge(x = typeIII_fixed_ames2nam_local[[2]], y = intersect_dong_ranges,
                             by.x = "rn2", by.y = "combined", all.x = TRUE)


typeIII_fixed_nam2nam_local[[2]]$rn2 <- gsub("_PC[0-9]", "", typeIII_fixed_nam2nam_local[[2]]$rn)
fixed_nam2nam_dong <- merge(x = typeIII_fixed_nam2nam_local[[2]], y = intersect_dong_ranges,
                             by.x = "rn2", by.y = "combined", all.x = TRUE)


# Gene window
gene_ames2nam_dong <- merge(x = typeIII_geneWindow_ames2nam_local[[2]], y = unique_gene,
                            by.x = "rn", by.y = "rn", all.x = TRUE)

gene_nam2nam_dong <- merge(x = typeIII_geneWindow_nam2nam_local[[2]], y = unique_gene,
                            by.x = "rn", by.y = "rn", all.x = TRUE)


# ---------------------------
# Plot out Type  III results
# ---------------------------

# Fixed size windows
ggplot(fixed_ames2nam_dong, aes(x = window.x, y = `Sum of Sq`, fill = PC)) +
  geom_bar(stat = 'identity') +
  geom_point(data = fixed_ames2nam_dong, aes(x = as.numeric(as.character(window.y)), y = pos), size = 1) +
  facet_grid(. ~ chr.x, scales = "free") +
  theme(legend.position="bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, 300)) +
  ggtitle("Type III SS from y(GDD_DTS) ~ 3 gPCs + 660 lPCs \n ames to nam applied PCs - 1 Mbp fixed window")

# Gene windows
ggplot(gene_ames2nam_dong, aes(x = window.x, y = `Sum of Sq`, fill = PC.x)) +
  geom_bar(stat = 'identity') +
  geom_point(data = gene_ames2nam_dong, aes(x = window.y, y = pos), size = 1) +
  facet_grid(. ~ chr, scales = "free") +
  theme(legend.position="bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, 300)) +
  ggtitle("Type III SS from y(GDD_DTS) ~ 3 gPCs + 660 lPCs \n ames to nam applied PCs - 181 genes/window")

# Fixed size windows
ggplot(fixed_nam2nam_dong, aes(x = window.x, y = `Sum of Sq`, fill = PC)) +
  geom_bar(stat = 'identity') +
  geom_point(data = fixed_nam2nam_dong, aes(x = as.numeric(as.character(window.y)), y = pos), size = 1) +
  facet_grid(. ~ chr.x, scales = "free") +
  theme(legend.position="bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, 300)) +
  ggtitle("Type III SS from y(GDD_DTS) ~ 3 gPCs + 660 lPCs \n nam only PCs - 1 Mbp fixed window")

# Gene windows
ggplot(gene_nam2nam_dong, aes(x = window.x, y = `Sum of Sq`, fill = PC.x)) +
  geom_bar(stat = 'identity') +
  geom_point(data = gene_ames2nam_dong, aes(x = window.y, y = pos), size = 1) +
  facet_grid(. ~ chr, scales = "free") +
  theme(legend.position="bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, 300)) +
  ggtitle("Type III SS from y(GDD_DTS) ~ 3 gPCs + 660 lPCs \n nam only PCs - 181 genes/window")


# ----------------------------------------
# Model R^2 values for numerous phenotypes
#    for each of the four PC methods
# ----------------------------------------

# Phenotype files
r_vals <- data.frame()
get_r2_models <- function(rowNames, globalPCs, localPCs, method_name){
  # Filter down phenotypes to only flowering time
  all_NAM_phenos <- all_NAM_phenos[,c(2,8:10,17:21,36:50,106:108,461:485)]
  
  # Combine global and local PCs, format taxa names
  nam_global_local <- cbind(rowNames, 
                            data.frame(globalPCs), 
                            data.frame(localPCs))
  
  # Change first column name
  colnames(nam_global_local)[1] <- "taxa"
  
  # Removes sequence identifiers from rownames
  nam_global_local[,1] <- gsub(":[0-9]{9}", "", nam_global_local[,1])
  
  # Merge dataframes together
  merged_data <- merge(x = all_NAM_phenos, y = nam_global_local,
                                 by.x = "Geno_Code", by.y = "taxa")
  
  # Separate phenotypes and PCs
  pcs <- merged_data[,53:715]
  phenos <- merged_data[,2:52]
  
  # Run model for all phenotypes
  for (i in seq(1:ncol(phenos))){
    model <- lm(phenos[,i] ~ ., pcs)
    tmp <- data.frame(colnames(phenos)[i], summary(model)$r.squared, summary(model)$adj.r.squared)
    tmp <- cbind(tmp, method_name)
    r_vals <- rbind(r_vals, tmp)
  }
  colnames(r_vals) <- c("phenotype", "r_squared", "adj_r_squared", "method")
  return(r_vals)
}

# Fixed window methods
fixed_ames2nam_r2 <- get_r2_models(rowNames = rownames(X_all$X_all), ames2nam_gPCs, fixed_ames2nam_local, "fixed_ames2nam")
fixed_nam2nam_r2 <- get_r2_models(rowNames = rownames(X_all_nam$X_all), nam2nam_gPCs, fixed_nam2nam_local, "fixed_nam2nam")

# gene window methods
geneWindow_ames2nam_r2 <- get_r2_models(rowNames = rownames(X_all$X_all), ames2nam_gPCs, geneWindow_ames2nam_local, "geneWindow_ames2nam")
geneWindow_nam2nam_r2 <- get_r2_models(rowNames = rownames(X_all_nam$X_all), nam2nam_gPCs, geneWindow_nam2nam_local, "geneWindow_nam2nam")

# COMBINE ALL
all_r2_across_methods <- rbind(fixed_ames2nam_r2, geneWindow_ames2nam_r2,
                               fixed_nam2nam_r2, geneWindow_nam2nam_r2)

# Plot put results
ggplot(all_r2_across_methods, aes(x=method, y=r_squared, fill=method)) + 
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  labs(title="R^2 from y(51 NAM phenotypes) ~ 3 gPCs + 660 lPCs(various methods)",
       x="Method", y = "R^2") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position="none")

ggplot(all_r2_across_methods, aes(x=method, y=adj_r_squared, fill=method)) + 
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  labs(title="Adjusted R^2 from y(51 NAM phenotypes) ~ 3 gPCs + 660 lPCs(various methods)",
       x="Method", y = "Adjusted R^2") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position="none")


# ---------------------------------------
# Look at Type III SS within Dong windows 
# vs. not in Dong windows for gene-windows
# ---------------------------------------
library(tidyverse)
# Function to subset data and calculate model comparison metrics
holder <- data.frame()
model_met <- function(dong_window_df, annotation){
  
  # Only take PC windows that contain a flowering time gene in v4
  windows_with_dongFT <- dong_window_df[!is.na(dong_window_df$v4_assoc_gene_model),]
  
  # Get all other windows with no flowering time gene in v4
  windows_without_dongFT <- dong_window_df[!which(dong_window_df$rn %in% windows_with_dongFT$rn),]
  
  # Compile differences in type III SS
  met1 <- sum(windows_with_dongFT$`Sum of Sq`)
  met2 <- sum(windows_with_dongFT$`Sum of Sq`)/nrow(windows_with_dongFT)
  met3 <- sum(windows_without_dongFT$`Sum of Sq`)/nrow(windows_without_dongFT)
  met4 <- (sum(windows_with_dongFT$`Sum of Sq`)/nrow(windows_with_dongFT))/(sum(windows_without_dongFT$`Sum of Sq`)/nrow(windows_without_dongFT))

  # Combine models and add helpful names
  holder <- rbind(holder, c(data.frame(annotation), met1, met2, met3, met4))
  colnames(holder) <- c("model_ID", 
                          "sum FT windows",
                          "sum within FT windows / number of FT windows",
                          "sum in non-FT windows/number of non-FT windows",
                          "sum within FT windows/sum non-FT windows")
  
  return(holder)
  
}


#  Use function to compare all methods
compare_ft_models <- rbind(model_met(fixed_ames2nam_dong, "fixed window - ames2nam"),
                           model_met(gene_ames2nam_dong, "gene window - ames2nam"),
                           model_met(fixed_nam2nam_dong, "fixed window - nam2nam"),
                           model_met(gene_nam2nam_dong, "gene window - nam2nam"))

print(compare_ft_models)


# -----------------------------
# Export data in tassel format
# -----------------------------

# Exporting 3 different things, steps
# - Subset phenotypes, remove missing NAM taxa
# - Export phenotypes in tassel format
# - Export global PCs in tassel format
# - Export global and local PCs in tassel format

# Load in NAM data
all_NAM_phenos <- read.csv("~/Box Sync/Cornell_PhD/labProjects/hackathons/2019-10-21_hackathon/all_NAM_phenos.txt", sep="")

# Filter down phenotypes to only 
# fumarate,malate,chlorophyll A,chlorophyll b, cobdiameter, southern leaf blight, flowering time
all_NAM_phenos <- all_NAM_phenos[,c(2,40,39,37,38,17,105,8:10)]

# Combine global and local PCs, format taxa names
nam_global_local <- cbind(rownames(X_all_nam$X_all), 
                          data.frame(nam2nam_gPCs), 
                          data.frame(geneWindow_nam2nam_local))

# Format PCs for merging
nam_global_local$copied_taxa <- gsub(":[0-9]{9}", "", nam_global_local[,1])
nam_global_local$copied_taxa <- gsub("_", "", nam_global_local$copied_taxa)

# Merge dataframes together
merged_pcs_phenos <- merge(x = all_NAM_phenos, y = nam_global_local,
                           by.x = "Geno_Code", by.y = "copied_taxa")

# Rearrage columns
merged_pcs_phenos <- merged_pcs_phenos[,c(11,2:10,12:ncol(merged_pcs_phenos))]

# Fix column  names
colnames(merged_pcs_phenos)[1] <- "taxa"

# Remove missing values
merged_pcs_phenos <- na.omit(merged_pcs_phenos)

#gemma --> stop here if you want to export in gemma format

# Make column names the first row *tassel format*
merged_pcs_phenos[,1] <- as.character(merged_pcs_phenos[,1])
merged_pcs_phenos <- rbind(colnames(merged_pcs_phenos), merged_pcs_phenos)

# Create covariate names
tassel_covar_and_data_names <- c("taxa", rep("data", ncol(all_NAM_phenos)-1), rep("covariate", ncol(nam_global_local)-1))

# Add in covariate, data, taxa IDs
merged_pcs_phenos <- rbind(tassel_covar_and_data_names, merged_pcs_phenos)


# Subset data based on tests

# Phenotypes only
phenos_only <- merged_pcs_phenos[,1:10]

# Global PCs only
global_only <- merged_pcs_phenos[,c(1,11:13)]

# Global and local PCs
global_local <- merged_pcs_phenos[c(1,11:ncol(merged_pcs_phenos))]

# Export data
setwd("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/data")
write.table(phenos_only, file = "range_complexity_phenos_wallace.txt",
            row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
write.table(global_only, file = "nam2nam_3gPCs.txt",
            row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
write.table(global_local, file = "nam2nam_3gPCs_660lPCs.txt",
            row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)


# Export in gemma format (run everything until #gemma indicator)
setwd("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/data")
write.table(merged_pcs_phenos[,2:10], "gemma_phenos_range_complexity.txt", 
            col.names = F, row.names = F, sep = " ", quote = F)
write.table(merged_pcs_phenos[,11:13], "gemma_covars_range_complexity.txt", 
            col.names = F, row.names = F, sep = " ", quote = F)

# Z006E0118:250022587 is duplicated in list

