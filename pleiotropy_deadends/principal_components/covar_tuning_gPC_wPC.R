# ---------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2019-12-19 
#
# Description 
#   - Parameter tuning for GWAS models
#   - Trying many models
      # method (ames2nam, nam2nam)
      # number of global PCs (3 or 26?)
      # number of local PCs (1-10?)
      # window size/window type
      #   - fixed window (500 kb - 1 Mbp)
      #   - gene window (100-180 in 10 gene intervals, currently at 180 genes/window)
      #   - genetic distance (200-400 windows, look at file to see how many cM this corresponds to)
      # Run all in Fast Association on same set of phenotypes
      # Benchmark to gPC+kinship model from Gemma (use SNP relate to go from vcf to Plink format)
# ---------------------------------------------------------------

# Load in source scripts
source('~/git_projects/haplotype_v_snp_gwas/src/R/pca_testing/ames2any_matrix.R')
source('~/git_projects/haplotype_v_snp_gwas/src/R/pca_testing/local_window_funs.R')

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
ames.vcf <- "~/Box Sync/Cornell_PhD/labProjects/hackathons/2019-10-21_hackathon/ZeaGBSv27_20171204_AGPv4_Ames_filteredGenotypes_beagleImputation_Filter_intersectSitesWithNAM_Filter.vcf"
ames_vcf <- snpgdsVCF2GDS(ames.vcf, "ames.gds", method="copy.num.of.ref", snpfirstdim=TRUE, ignore.chr.prefix = "S")
snpgdsSummary("ames.gds")
genofile_ames <- snpgdsOpen(ames_vcf)

# Load in NAM SNPs that are shared (intersect) with the Ames SNPs (did intersection in TASSEL)
# REPLACE ME WITH IBM INCLUDED SNPS/TAXA
nam.vcf <- "~/Box Sync/Cornell_PhD/labProjects/hap_gwa/data/genos/intersectJoin_NAMbyAmes_shared40k_allNAM_IBM_beagleImputation.vcf.gz"
nam_vcf <- snpgdsVCF2GDS(nam.vcf, "nam.gds", method="copy.num.of.ref", snpfirstdim=TRUE, ignore.chr.prefix = "S")
snpgdsSummary("nam.gds")
genofile_nam <- snpgdsOpen(nam_vcf)

# Load in nam founder data
nam_founder.vcf <- "~/Box Sync/Cornell_PhD/labProjects/hap_gwa/data/genos/intersectJoin_NAMbyAmes_shared_beagleImputed_NAMfounders.vcf.gz"
nam_founder_vcf <- snpgdsVCF2GDS(nam_founder.vcf, "nam_founder.gds", method="copy.num.of.ref", snpfirstdim=TRUE, ignore.chr.prefix = "S")
snpgdsSummary("nam_founder.gds")
genofile_nam_founder <- snpgdsOpen(nam_founder_vcf)


# ------------------------------------------
#     Calculate ames2nam global PCs
#     uses ames2any_matrix functions
# ------------------------------------------

# Combine genofile objects
X_allNAM <- combine_genofile_matrix(genofile_ames = genofile_ames, genofile_any = genofile_nam)

# ------
# 3 gPCs
# ------ 

# Calculate Ames global PCs, get loadings
global_ames_BV_3gPCs <- get_ames_BV(X_ames = genofile_ames, Q = NULL, num_PCs = 3)

# Transfer loadings and coefficients over to all X (Ames, NAM) to get adjusted global PCs
# Get Q matrix (population structure)
ames2nam_gPCs_3gPCs_allNAM <- ames2any(X_all = X_allNAM$X_all,
                          Q = NULL, B = global_ames_BV_3gPCs$B, V = global_ames_BV_3gPCs$V)
Q_all_3gPCs_allNAM <- cbind(1, ames2nam_gPCs_3gPCs_allNAM)

# ------
# 26 gPCs
# ------ 

# Calculate Ames global PCs, get loadings
global_ames_BV_26gPCs <- get_ames_BV(X_ames = genofile_ames, Q = NULL, num_PCs = 26)
ames2nam_gPCs_26gPCs_allNAM <- ames2any(X_all = X_allNAM$X_all,
                                 Q = NULL, B = global_ames_BV_26gPCs$B, V = global_ames_BV_26gPCs$V)
Q_all_26gPCs_allNAM <- cbind(1, ames2nam_gPCs_26gPCs_allNAM)


# ------------------------------------------
#    Calculate local PCs in fixed windows
# ------------------------------------------


# Use fixed widow function on ames to nam calculated PCs
# 3 gPCs
fixed_ames2nam_local_3allNAM <- calc_local_pcs_fixed_window(X_all = X_allNAM,
                                                    genofileID_ames = genofile_ames,
                                                    Q_all = Q_all_3gPCs_allNAM,
                                                    numLocalPcs = 3, windowSize = 1e7)
# 26 gPCs
fixed_ames2nam_local_26allNAM <- calc_local_pcs_fixed_window(X_all = X_allNAM,
                                                            genofileID_ames = genofile_ames,
                                                            Q_all = Q_all_26gPCs_allNAM,
                                                            numLocalPcs = 3, windowSize = 1e7)


# ------------------------------------------
#   Calculate local PCs in defined windows
#           formatting step
# ------------------------------------------

# Get all maize v4 genes from gff file 
maize_gff <- read.gff("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/data/maize_v4_annotations/Zea_mays.B73_RefGen_v4.45.gff3", 
                      na.strings = c(".", "?"), GFF3 = TRUE)

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
# 3 gPCs
geneWindow_ames2nam_local_3allNAM <- calc_local_pcs_gene_window(X_all = X_allNAM,
                                                        genofileID_ames = genofile_ames,
                                                        Q_all = Q_all_3gPCs_allNAM,
                                                        numLocalPcs = 3,
                                                        windowFile = window_file)

# 26 gPCs
geneWindow_ames2nam_local_26allNAM <- calc_local_pcs_gene_window(X_all = X_allNAM,
                                                        genofileID_ames = genofile_ames,
                                                        Q_all = Q_all_26gPCs_allNAM,
                                                        numLocalPcs = 3,
                                                        windowFile = window_file)


# -----------------------------
# Export files in tassel format
# -----------------------------

# Load in NAM data
# Filter down phenotypes to only 
# fumarate,malate,chlorophyll A,chlorophyll b, cobdiameter, southern leaf blight, flowering time
all_NAM_phenos <- read.csv("~/Box Sync/Cornell_PhD/labProjects/hackathons/2019-10-21_hackathon/all_NAM_phenos.txt", sep="")
all_NAM_phenos <- all_NAM_phenos[,c(2,4,40,39,37,38,17,105,8:10)]

# "Fix" IBM names
all_NAM_phenos$IBM_compat <- as.character(all_NAM_phenos$Geno_Code)
all_NAM_phenos[c(3203:3402),12] <- as.character(all_NAM_phenos[c(3203:3402),2])
all_NAM_phenos$IBM_compat <- gsub("MO", "M0", all_NAM_phenos$IBM_compat)

# Function to format data
tasselize <- function(main_matrix, global_matrix, window_matrix, filename){
  # Combine objects
  together <- cbind(data.frame(rownames(main_matrix)), data.frame(global_matrix), data.frame(window_matrix))

  # Format taxa IDs from PCs for merging
  together$copied_taxa <- gsub(":[0-9]{9}", "", together[,1])
  together$copied_taxa <- gsub("_", "", together$copied_taxa)
  print(str(together))

  # Merge dataframes together
  together <- merge(x = all_NAM_phenos, y = together, sort = F,
                             by.x = "IBM_compat", by.y = "copied_taxa")
  print("post merge")
  print(str(together))
  # Sort columns (3 or 26 gPCs)
  together <- together[,c(13,4:12,14:ncol(together))]

  # Change name of first column from taxa ID to trait
  colnames(together)[1] <- "<Trait>"
  
  # Remove missing and duplicated values
  together <- na.omit(together)
  together <- together[!duplicated(together),]
  
  # To export phenotypes as well
  # return(together)
  
  # Make column names the first row
  together[,1] <- as.character(together[,1])
  header <- c("<Covariate>", rep("", ncol(together)-1))
  together <- rbind(header, colnames(together), together)

  # Only get PCs
  together <- together[,c(1,11:ncol(together))]

  # Export
  write.table(together, file = filename, row.names = FALSE, 
              col.names = FALSE, sep = "\t", quote = FALSE)

}

# fixed window options
tasselize(X_allNAM$X_all, ames2nam_gPCs_3gPCs_allNAM, fixed_ames2nam_local_3allNAM, 
          filename = "ames2nam_3gPCs_allNAM_fixedWindow.txt")
tasselize(X_allNAM$X_all, ames2nam_gPCs_26gPCs_allNAM, fixed_ames2nam_local_26allNAM, 
          filename = "ames2nam_26gPCs_allNAM_fixedWindow.txt")

# Gene window options
tasselize(X_allNAM$X_all, ames2nam_gPCs_3gPCs_allNAM, geneWindow_ames2nam_local_3allNAM, 
          filename = "ames2nam_3gPCs_allNAM_geneWindow.txt")
tasselize(X_allNAM$X_all, ames2nam_gPCs_26gPCs_allNAM, geneWindow_ames2nam_local_26allNAM, 
          filename = "ames2nam_26gPCs_allNAM_geneWindow.txt")


# ----------------------------------------
#      Calculate nam2nam global PCs
#         (no Ames adjustment)
# ----------------------------------------

nam2nam_test_methods <- function(start_matrix, end_matrix, num_gPCs = 3, num_wPCs = 3, out_filename){
  
  # Calculate NAM global B and V matrix
  global_nam_BV <- get_ames_BV(X_ames = start_matrix, Q = NULL, num_PCs = num_gPCs)
  
  # Make NAM X matrix and do formatting
  X_all_nam_geno <- snpgdsGetGeno(end_matrix, snpfirstdim = FALSE, with.id = TRUE, verbose = FALSE)
  X_all_nam <- 2-X_all_nam_geno$genotype
  colnames(X_all_nam) <- X_all_nam_geno$snp.id
  rownames(X_all_nam) <- c(X_all_nam_geno$sample.id)
  X_all_nam <- list(X_all = X_all_nam, snp_info = snpgdsSNPList(end_matrix))
  
  # Calculate PCs
  nam2nam_gPCs <- ames2any(X_all = X_all_nam$X_all,
                           Q = NULL, B = global_nam_BV$B, V = global_nam_BV$V)
  
  # Make Q matrix (population structure)
  Q_all_nam <- cbind(1, nam2nam_gPCs)
  print("Finshed calculating and formatting global PCs")
  
  # Use fixed window function
  fixed_nam2nam_local <- calc_local_pcs_fixed_window(X_all = X_all_nam,
                                                     genofileID_ames = start_matrix,
                                                     Q_all = Q_all_nam,
                                                     numLocalPcs = num_wPCs,
                                                     windowSize = 1e7)
  print("Finshed calculating fixed window PCs")
  
  # Use gene window function
  geneWindow_nam2nam_local <- calc_local_pcs_gene_window(X_all = X_all_nam,
                                                         genofileID_ames = start_matrix,
                                                         Q_all = Q_all_nam,
                                                         numLocalPcs = num_wPCs,
                                                         windowFile = window_file)
  print("Finshed calculating gene window PCs")
  
  # Tasselize and export
  filenameFixed <- paste(out_filename, "_fixedWindow.txt", sep = "")
  filenameGene <- paste(out_filename, "_geneWindow.txt", sep = "")
  
  tasselize(main_matrix = X_all_nam$X_all, global_matrix = nam2nam_gPCs,
            window_matrix = fixed_nam2nam_local, filename = filenameFixed)
  
  tasselize(main_matrix = X_all_nam$X_all, global_matrix = nam2nam_gPCs,
            window_matrix = geneWindow_nam2nam_local, filename = filenameGene)
  print("Finshed exporting both files")
  
  return(list(global = nam2nam_gPCs, fixed = fixed_nam2nam_local, gene = geneWindow_nam2nam_local))
}

# nam 2 nam all-taxa methods
allnam2allnam_3gPCs <- nam2nam_test_methods(genofile_nam, genofile_nam, num_gPCs = 3, out_filename = "allnam2allnam_3gPCs")
allnam2allnam_26gPCs <- nam2nam_test_methods(genofile_nam, genofile_nam, num_gPCs = 26, out_filename = "allnam2allnam_26gPCs")

# Nam 2 nam founders methods
namfounder2allnam_3gPCs <- nam2nam_test_methods(genofile_nam_founder, genofile_nam, num_gPCs = 3, out_filename = "namfounder2allnam_3gPCs")
namfounder2allnam_26gPCs <- nam2nam_test_methods(genofile_nam_founder, genofile_nam, num_gPCs = 26, out_filename = "namfounder2allnam_26gPCs")



