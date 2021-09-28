# ---------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2020-03-30
#
# Description 
#   - Script to run rTASSEL on window PCs, dropping out windows
#   - in analyzed SNP's region
# ---------------------------------------------------------------

# Set working directory
setwd("/workdir/mbb262/mer_debugging")

# Install package
if (!require("devtools")) install.packages("devtools")
devtools::install_bitbucket(
  repo = "bucklerlab/rtassel",
  ref = "rTAS-006-genome-filtration",
  build_vignettes = F
)

# Setting memory
options(java.parameters = c("-Xmx200g"))

# Load package
library(rTASSEL)

# Start logging
rTASSEL::startLogger(fullPath = NULL, fileName = NULL)


# -------------
#   Genotypes
# -------------

# Load chromosome 8 and 10 vcf data
chr8_tasGenoVCF <- rTASSEL::readGenotypeTableFromPath(
  path = "/home/mbb262/tassel_gwa/data/genos/combined_ibm_nam_all_imputed/nam_ibm_imputed_intersectMAFFilter_chr8.vcf.gz",
  keepDepth = FALSE)
chr8_tasGenoVCF

chr10_tasGenoVCF <- rTASSEL::readGenotypeTableFromPath(
  path = "/home/mbb262/tassel_gwa/data/genos/combined_ibm_nam_all_imputed/nam_ibm_imputed_intersectMAFFilter_chr10.vcf.gz",
  keepDepth = TRUE)
chr10_tasGenoVCF


# --------------
#   Phenotypes
# --------------

# Read from phenotype path
# Load into rTASSEL `TasselGenotypePhenotype` object
wallace_phenos <- rTASSEL::readPhenotypeFromPath(path = "/home/mbb262/tassel_gwa/data/phenos/2020_03_24_variablewPCs/wallace_phenos_2020_03_16.txt")

# Inspect object
# Traits: Taxa cob_diameter_raw_Hung2012 Days_To_Silk_BLUP_Sum0607_Buckler2009 
wallace_phenos

# Turn into R dataframe
wallace_phenosExportPhenoDF <- rTASSEL::getPhenotypeDF(tasObj = wallace_phenos)
wallace_phenosExportPhenoDF


# --------------
#   Covariates
# --------------

# Load in global PC path
# Load into rTASSEL `TasselGenotypePhenotype` object
gPC_Path  <- "/home/mbb262/tassel_gwa/data/phenos/2020_03_24_variablewPCs/ames2nam_3gPCs_allNAM_gPCs_only.txt"
gPCs <- rTASSEL::readPhenotypeFromPath(path = gPC_Path)
gPCs

# Turn into R dataframe
gPCExportPhenoDF <- rTASSEL::getPhenotypeDF(tasObj = gPCs)
gPCExportPhenoDF

# Load in two rounds of window PCs
main_wPCs <- read.csv("/home/mbb262/tassel_gwa/data/phenos/2020_03_24_variablewPCs/mainGeneWindow360_ames2nam_local_3allNAM_varaibleNumPCs.csv", header = TRUE)
# midway_wPCs <- read.csv("/home/mbb262/tassel_gwa/data/phenos/2020_03_24_variablewPCs/midwayGeneWindow360_ames2nam_local_3allNAM_varaibleNumPCs.csv", header = TRUE)
midway_wPCs <- read.csv("/home/mbb262/tassel_gwa/data/phenos/2020_03_24_variablewPCs/midwayGeneWindow360_ames2nam_local_3allNAM_varaibleNumPCs_2.csv", header = TRUE)


# ---------------------
# Merge files together
# ---------------------

# Merge phenotypes and global PCs
pheno_gpcs <- merge(x = wallace_phenosExportPhenoDF, y = gPCExportPhenoDF)


# ---------------
#  Other Stuff
# ---------------

# Load in window objects
window_file_midway_360 <- read.csv("/home/mbb262/tassel_gwa/data/phenos/2020_03_24_variablewPCs/window_file_midway_360.csv")


# ------------------
# Debugging
# -----------------
library(dplyr)
window_file_df <- window_file_midway_360 %>% filter(chrom == 8, windowLocation == "mainWindow") %>% arrange(chrom, start)
window <- 4
window_file_df[window,]

tasGenoFilt <- rTASSEL::filterGenotypeTableSites(tasObj = chr8_tasGenoVCF,
                                                 siteRangeFilterType = "position",
                                                 startChr = 8, endChr = 8,
                                                 startPos = window_file_df$start[window],
                                                 endPos = window_file_df$stop[window])

# Subtract out window PCs within this region
# Command subtacts out columns with grep matches of window name, needed to add "_" to avoid partial matches
window_pcs_df <- main_wPCs
print(ncol(window_pcs_df))
subset_wPCs <- window_pcs_df[ , -grep(paste(window_file_df$windowID[window], "_", sep = ""), colnames(window_pcs_df))]
print(ncol(subset_wPCs))
# Print reassuring messages
message(paste("We want to drop this column: ", window_file_df$windowID[window]))
message("columns dropped in this window")
print(colnames(window_pcs_df[ , grep(paste("^", window_file_df$windowID[window], sep = ""), colnames(window_pcs_df))]))

# Intersect genotypes, phenotypes, and PCs
pheno_gpcs_windowPCs <- merge(x = pheno_gpcs, y = subset_wPCs, by.x = "Taxa", by.y = "X")

tasselize <- rbind(c("<Phenotype>", rep("", ncol(pheno_gpcs_windowPCs)-1)), 
                   c("taxa", "data", "data", rep("covariate", ncol(pheno_gpcs_windowPCs)-3)), 
                   colnames(pheno_gpcs_windowPCs), 
                   pheno_gpcs_windowPCs)
write.table(tasselize, "chr8_window4_pcs.txt", col.names = F, row.names = F, quote = F, sep = "\t")


# Tasselize merged phenotypes + global PCs + window PCs
tasPhenoDF <- rTASSEL::readPhenotypeFromDataFrame(
  phenotypeDF = pheno_gpcs_windowPCs,
  taxaID = "Taxa",
  attributeTypes = c("data", "data", rep("covariate", ncol(pheno_gpcs_windowPCs)-3)))

# Join genotypes with phenotypes + g PCs + wPCs
tasPhenoDF <- rTASSEL::readGenotypePhenotype(
  genoPathOrObj = tasGenoFilt,
  phenoPathDFOrObj = tasPhenoDF)
print(tasPhenoDF)



# Run Fast Associaion
tasFAST <- rTASSEL::assocModelFitter(
  tasObj = tasPhenoDF,
  formula = Days_To_Silk_BLUP_Sum0607_Buckler2009 ~ .,
  fitMarkers = TRUE,
  kinship = NULL,
  fastAssociation = TRUE,
  maxP = 0.0001,
  maxThreads = 100)
print(tasFAST)



# ------------------------
# Start Analysis Pipeline
# ------------------------

# Function subsets a TASSEL genotype table, subsets it, and subets window PCs
#   then runs fast association
# Function inputs a Tassel genotype object 
library(dplyr)
run_fast_assoc_dropout <- function(tassel_genotype, window_file_df, window_pcs_df, tas_phenos_gPCs, chr, windowType = c("mainWindow", "midwayWindow")){
  
  # Subset windowfile results on chromosome
  window_file_df <- window_file_df %>% filter(chrom == chr, windowLocation == windowType) %>% arrange(chrom, start)
  print(window_file_df)
  # Collect results outside of loop
  fa_results <- c()
  
  # Start the loop
  for (window in seq_len(nrow(window_file_df))){
    
    # Print status message
    message(paste("I am on window:", window, " out of ", nrow(window_file_df)))
    
    # Subset genotypes on the start and stop position of single item in window file
    tasGenoFilt <- rTASSEL::filterGenotypeTableSites(tasObj = tassel_genotype,
                                                     siteRangeFilterType = "position",
                                                     startChr = chr, endChr = chr,
                                                     startPos = window_file_df$start[window],
                                                     endPos = window_file_df$stop[window])
    
    # Subtract out window PCs within this region
    # Command subtacts out columns with grep matches of window name, needed to add "_" to avoid partial matches
    print(ncol(window_pcs_df))
    subset_wPCs <- window_pcs_df[ , -grep(paste(window_file_df$windowID[window], "_", sep = ""), colnames(window_pcs_df))]
    print(ncol(subset_wPCs))
    # Print reassuring messages
    message(paste("We want to drop this column: ", window_file_df$windowID[window]))
    message("columns dropped in this window")
    print(colnames( window_pcs_df[ , grep(paste("^", window_file_df$windowID[window], sep = ""), colnames(window_pcs_df))]))
    
    # Intersect genotypes, phenotypes, and PCs
    pheno_gpcs_windowPCs <- merge(x = tas_phenos_gPCs, y = subset_wPCs, by.x = "Taxa", by.y = "X")
    
    # Tasselize merged phenotypes + global PCs + window PCs
    tasPhenoDF <- rTASSEL::readPhenotypeFromDataFrame(
      phenotypeDF = pheno_gpcs_windowPCs,
      taxaID = "Taxa",
      attributeTypes = c("data", "data", rep("covariate", ncol(pheno_gpcs_windowPCs)-3)))
    
    # Join genotypes with phenotypes + g PCs + wPCs
    tasPhenoDF <- rTASSEL::readGenotypePhenotype(
      genoPathOrObj = tasGenoFilt,
      phenoPathDFOrObj = tasPhenoDF)
    print(tasPhenoDF)
    
    # Debug run fast association
    assoc_bruh <- jRC$fastAssociation(
      tasPhenoDF@jTasselObj,
      rJava::.jnew("java/lang/Double", 1),
      rJava::.jnew("java/lang/Integer", "35")
    )
    fast_tab_object <- assoc_bruh$get("FastAssociation")
    fast_tab_object$getRowCount()
    fast_tab_object$getElementCount()
    fast_tab_object$getValueAt(rJava::.jlong(0), as.integer(0))
    brandon <- debugConvert(fast_tab_object)
    print(brandon)
    print(head(brandon))
    fa_results <- rbind(fa_results, brandon)
    # Run Fast Associaion
    # tasFAST <- rTASSEL::assocModelFitter(
    #   tasObj = tasPhenoDF,
    #   formula = Days_To_Silk_BLUP_Sum0607_Buckler2009 ~ .,
    #   fitMarkers = TRUE,
    #   kinship = NULL,
    #   fastAssociation = TRUE,
    #   maxP = 0.0001,
    #   maxThreads = 100)
    # print(tasFAST)
    
    # Combine this window with the prior windows
    # names <- data.frame(rep(window_file_df$windowLocation[window], nrow(tasFAST$FastAssociation)))
    # winName <- data.frame(rep(window_file_df$windowID[window], nrow(tasFAST$FastAssociation)))
    # temp <- cbind(tasFAST$FastAssociation, names, winName)
    # colnames(temp)[8] <- "windowLocation"
    # colnames(temp)[9] <- "windowID"
    # fa_results <- rbind(fa_results, temp)
  }
  # Return results
  return(fa_results)
}

# Use function
chr8_main <- run_fast_assoc_dropout(tassel_genotype = chr8_tasGenoVCF,
                                    window_file_df = window_file_midway_360,
                                    window_pcs_df = main_wPCs,
                                    tas_phenos_gPCs = pheno_gpcs,
                                    chr = 8, windowType = "mainWindow")
write.csv(chr8_main, "chr8_main_faresults.csv", row.names = F)

chr10_main <- run_fast_assoc_dropout(tassel_genotype = chr10_tasGenoVCF, 
                                     window_file_df = window_file_midway_360,
                                     window_pcs_df = main_wPCs,
                                     tas_phenos_gPCs = pheno_gpcs,
                                     chr = 10, windowType = "mainWindow")
write.csv(chr10_main, "chr10_main_faresults.csv", row.names = F)
chr10_main <- read.csv("chr10_main_faresults.csv", header = TRUE)

# Forgive me jesus I hardcoded names in the midway_wPC files for chrom 8 and 10 --> will change later
chr10_midway <- run_fast_assoc_dropout(tassel_genotype = chr10_tasGenoVCF,
                                       window_file_df = window_file_midway_360,
                                       window_pcs_df = midway_wPCs,
                                       tas_phenos_gPCs = pheno_gpcs,
                                       chr = 10, windowType = "midwayWindow")
write.csv(chr10_midway, "chr10_midway_faresults.csv", row.names = F)
chr10_midway <- read.csv("chr10_midway_faresults.csv", header = TRUE)

# Made a mistake when naming window PCs, manually change them
chr8_midway <- run_fast_assoc_dropout(tassel_genotype = chr8_tasGenoVCF,
                                      window_file_df = window_file_midway_360,
                                      window_pcs_df = midway_wPCs,
                                      tas_phenos_gPCs = pheno_gpcs,
                                      chr = 8, windowType = "midwayWindow")
write.csv(chr8_midway, "chr8_midway_faresults.csv", row.names = F)

# Merge files
# Merge datasets by marker
main_midway_chr8 <- merge(chr8_main, chr8_midway, by = "Marker", suffixes = c("_main", "_midway"), all = TRUE)
main_midway_chr10 <- merge(chr10_main, chr10_midway, by = "Marker", suffixes = c("_main", "_midway"), all = TRUE)

# Find min p-value row (out of two columns) and save a new column with lowest p-value
main_midway_chr8 <- transform(main_midway_chr8, p = pmin(p_main, p_midway, na.rm = TRUE))
main_midway_chr10 <- transform(main_midway_chr10, p = pmin(p_main, p_midway, na.rm = TRUE))

# Change column names to be easier to plot
colnames(main_midway_chr8)[2:4] <- c("Trait", "Chr", "Pos")
colnames(main_midway_chr10)[2:4] <- c("Trait", "Chr", "Pos")

# Load in old data with all PCs applied
main_m10_chr8 <- read.table("~/tassel_gwa/results/2020_03_28_mainAndMidwayWindows/model10_chrom8_ames2nam_3gPCs_main360geneWindow_FAresults_2020_03_28.txt", header = TRUE)
main_m10_chr10 <- read.table("~/tassel_gwa/results/2020_03_28_mainAndMidwayWindows/model10_chrom10_ames2nam_3gPCs_main360geneWindow_FAresults_2020_03_28.txt", header = TRUE)
midway_m10_chr8 <- read.table("~/tassel_gwa/results/2020_03_28_mainAndMidwayWindows/model10_chrom8_ames2nam_3gPCs_midway360geneWindow_FAresults_2020_03_28.txt", header = TRUE)
midway_m10_chr10 <- read.table("~/tassel_gwa/results/2020_03_28_mainAndMidwayWindows/model10_chrom10_ames2nam_3gPCs_midway360geneWindow_FAresults_2020_03_28.txt", header = TRUE)

# Subset out only DTS
main_m10_chr8 <- main_m10_chr8 %>% filter(Trait == "Days_To_Silk_BLUP_Sum0607_Buckler2009")
main_m10_chr10 <- main_m10_chr10 %>% filter(Trait == "Days_To_Silk_BLUP_Sum0607_Buckler2009")
midway_m10_chr8 <- midway_m10_chr8 %>% filter(Trait == "Days_To_Silk_BLUP_Sum0607_Buckler2009")
midway_m10_chr10 <- midway_m10_chr10 %>% filter(Trait == "Days_To_Silk_BLUP_Sum0607_Buckler2009")

# Merge datasets by factorized position
all_main_midway_chr8 <- merge(main_m10_chr8, midway_m10_chr8, by = "Marker", suffixes = c("_main", "_midway"))
all_main_midway_chr10 <- merge(main_m10_chr10, midway_m10_chr10, by = "Marker", suffixes = c("_main", "_midway"))

# Find min p-value row (out of two columns) and save a new column with lowest p-value
all_main_midway_chr8$p <- apply(all_main_midway_chr8[,c(7,13)],1,min,na.rm=TRUE)
all_main_midway_chr10$p <- apply(all_main_midway_chr10[,c(7,13)],1,min,na.rm=TRUE)

# Rename columns
colnames(all_main_midway_chr8)[2:4] <- c("Trait", "Chr", "Pos")
colnames(all_main_midway_chr10)[2:4] <- c("Trait", "Chr", "Pos")

# Rank correlation
temp <- na.omit(main_midway_chr10[,c(7,15)])
cor.test( ~ p_main + p_midway,
          data=temp,
          method = "spearman",
          continuity = FALSE,
          conf.level = 0.95)
# Count how many are identical.
temp <- summarise(group_by(temp,-log10(p_main), -log10(p_midway)),length(p_midway))
temp <- group_size(group_by(temp,p_main, p_midway))

# -------------------------
#     Plot what I have
# -------------------------

source('~/git_projects/haplotype_v_snp_gwas/src/R/pca_testing/create_gene_windows.R')
a_chr8_drop <- manhattan_with_windows8(results = main_midway_chr8,
                                       upper_bound = 30,
                                       chr = 8,
                                       window_file_df = window_file_midway_360,
                                       title ="Window dropout - Chr 8 \n DTS ~ 3gPCs + 218 overlapping windows * 1-8 PCs per window")

b_chr10_drop <- manhattan_with_windows10(results = main_midway_chr10,
                                         upper_bound = 90,
                                         chr = 10,
                                         window_file_df = window_file_midway_360,
                                         title ="Window dropout - Chr 10 \n DTS ~ 3gPCs + 218 overlapping windows * 1-8 PCs per window")

b_chr8_all <- manhattan_with_windows8(results = all_main_midway_chr8,
                                      upper_bound = 30,
                                      chr = 8,
                                      window_file_df = window_file_midway_360,
                                      title ="All windows - Chr 8 \n DTS ~ 3gPCs + 218 overlapping windows * 1-8 PCs per window")
b_chr10_all <- manhattan_with_windows10(results = all_main_midway_chr10,
                                        upper_bound = 90,
                                        chr = 10,
                                        window_file_df = window_file_midway_360,
                                        title ="All windows - Chr 10 \n DTS ~ 3gPCs + 218 overlapping windows * 1-8 PCs per window")



library(patchwork)
(a_chr8_drop/b_chr8_all) | (b_chr10_all/b_chr10_drop)
ggsave("main_midway_window_dropout_and_all_comparison.png")


# ---------------------
# Rank correlation
# ---------------------

# Intersection of sites with p-values <0.0001 in two models
temp <- na.omit(main_midway_chr10[, c("p_main", "p_midway")])

cor.test( ~ p_main + p_midway,
          data=temp,
          method = "spearman",
          continuity = FALSE,
          conf.level = 0.95)
# Spearman's rank correlation rho
# 
# data:  p_main and p_midway
# S = 2.7552e+14, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#      rho 
# 0.758165 
plot(p_midway~p_main,
     data=temp,
     pch=16)


# Accessing elements directly through Java
# # rUn in hack-y way
jRC <- rJava::J("net/maizegenetics/plugindef/GenerateRCode")
debugConvert <- function(stringTab) {
  stringTab <- stringTab$toStringTabDelim()
  obj <- unlist(strsplit(stringTab, split = "\n"))
  obj <- strsplit(obj, split = "\t")
  obj <- t(simplify2array(obj))
  colnames(obj) <- as.character(unlist(obj[1, ]))
  obj <- obj[-1, ]
  # Check if object becomes named vector
  if (is.vector(obj)) {
    obj <- dplyr::bind_rows(obj)
    S4Vectors::DataFrame(obj)
  } else {
    S4Vectors::DataFrame(obj)
  }
}

assoc_bruh <- jRC$fastAssociation(
  tasPhenoDF@jTasselObj,
  rJava::.jnew("java/lang/Double", 1),
  rJava::.jnew("java/lang/Integer", "35")
)
fast_tab_object <- assoc_bruh$get("FastAssociation")
fast_tab_object$getRowCount()
fast_tab_object$getElementCount()
fast_tab_object$getValueAt(rJava::.jlong(0), as.integer(0))
brandon <- debugConvert(fast_tab_object)
