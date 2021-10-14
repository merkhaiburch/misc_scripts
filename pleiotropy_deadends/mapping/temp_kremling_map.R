
# Setting memory
options(java.parameters = c("-Xmx501g"))
options(java.parameters = c("-Xmx501g"))
options(java.parameters = c("-Xmx501g"))


# Load package
library(rTASSEL)
library(data.table)
library(magrittr)
library(dplyr)

# Start logging
rTASSEL::startLogger(fullPath = "/workdir/mbb262/results", 
                     fileName = NULL)

# Create global variables for mapping parameters
nThreads = 88
pValThresh = 1.0

# ------------------------------------------------------------------------------------
#                           Analysis pipeline Kremling
# ------------------------------------------------------------------------------------


# Set directory
setwd("/workdir/mbb262/results/kremling/")

# ----------------------------------
#   Load in pre-analyzed Covariates
# ----------------------------------

# Load in global PCs (Ames to Goodman method)
ames2goodman_gPCs_3gPCs <- read.csv("/home/mbb262/git_projects/pleiotropy/data/Q_all_3gPCs_allGoodman.csv", header = TRUE)
ames2goodman_gPCs_3gPCs <- ames2goodman_gPCs_3gPCs[,-c(5:6)]
colnames(ames2goodman_gPCs_3gPCs) <- c("Taxa", "PC1", "PC2", "PC3")

# Download in PEERs from irods
# icd /iplant/home/shared/commons_repo/curated/Kremling_Nature3RNASeq282_March2018/Expression_matrix/TASSEL_fmt_expression_w_covariates
# iget -r ./
# Remove first two lines with <Phenotype> from tassel header (I would load these into R with Rtassel but I keep getting an error)
# sed -i '1d' ./*.txt # ran twice because I'm not very clever
# sed -i '1d' ./*.txt

# Or, get PEERS from my blfs1 account --> do not have to remove first lines of file with sed commands
# mkdir /workdir/mbb262/peers_copy
# scp -r mbb262@cbsublfs1.tc.cornell.edu:/data1/users/mbb262/phenotypes/1_individual_datasets/kremling_2018_naturegenetics_expression/peers/*.txt /workdir/mbb262/peers_copy

# Load in peers
setwd("/workdir/mbb262/peers_copy/")
temp = list.files(pattern="*.txt")
myfiles = lapply(temp, 
                 data.table::fread, 
                 header=TRUE, 
                 sep = "\t", 
                 select = c(1, 7:31))

# Create shorter name for lists and then rename myfiles
short_names <- gsub("TASSEL_HEADER_df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames_", "", temp)# Merge the files together
short_names <- gsub("_sm_rand_and_box_coxed_w_MDSPCs_and_25PEERs_avgDups.txt", "", short_names)
names(myfiles) <- short_names

# Join with my PCs
myfiles <- lapply(myfiles, inner_join, x = ames2goodman_gPCs_3gPCs)


# ------------------------
# Format Kremling data
# ------------------------

# Get all file names
setwd("/workdir/mbb262/phenotypes/kremling/")
temp = list.files("/workdir/mbb262/phenotypes/kremling/", pattern="*.csv")

# Load data into r as sepearte dfs within a single list
formatted_kremling <- lapply(temp, 
                             data.table::fread, 
                             header=TRUE,
                             nThread = 80)

# Shorten names
short_names <- gsub("_kremling_formatted_v4_hapmapids.csv", "", temp)
names(formatted_kremling) <- short_names # rename each element in the list to the tissue

# Remove L3 from the list
formatted_kremling <- formatted_kremling[-4]

# Add on tissue name to gene names
for (i in seq_along(formatted_kremling)){
  colnames(formatted_kremling[[i]]) <- paste0(colnames(formatted_kremling[[i]]), "_", names(formatted_kremling[i]))
  colnames(formatted_kremling[[i]])[1] <- "Taxa"
}

# Check
lapply(formatted_kremling, function(x) x[1:5,1:5])


# -------------------------------------
# Join expression data with PCs + PEERs
# -------------------------------------

# Check order of both lists, must be the same
names(formatted_kremling)
names(myfiles)

# Innerjoin two lists by index and by taxa name
expression_pcs_peers <- purrr::map2(myfiles, formatted_kremling, inner_join, by = "Taxa")

# Subset out the three tissues we need
# LMAD <- expression_pcs_peers$LMAD
GShoot <- expression_pcs_peers$GShoot
# L3Base <- expression_pcs_peers$L3Base


# Tasselize merged phenotypes + PCs + PEERS
# LMAD
# tasPhenoDF_LMAD <- rTASSEL::readPhenotypeFromDataFrame(
#   phenotypeDF = LMAD,
#   taxaID = "Taxa",
#   attributeTypes = c(rep("covariate", 28), rep("data", ncol(LMAD)-29)))

# GShoot
tasPhenoDF_GShoot <- rTASSEL::readPhenotypeFromDataFrame(
  phenotypeDF = GShoot,
  taxaID = "Taxa",
  attributeTypes = c(rep("covariate", 28), rep("data", ncol(GShoot)-29)))

# # L3Base
# tasPhenoDF_L3Base <- rTASSEL::readPhenotypeFromDataFrame(
#   phenotypeDF = L3Base,
#   taxaID = "Taxa",
#   attributeTypes = c(rep("covariate", 28), rep("data", ncol(L3Base)-29)))



# ---------------------------
#       Map phenotypes
# RUN RTASSEL FOR EACH OF THE THREE ERRORED FILES
# ---------------------------

# Set directory for results
setwd("/workdir/mbb262/results/kremling/")
vcfpath_goodman <- "/workdir/mbb262/genotypes/goodman/"

# chrom 1, tissues GShoot and L3Base
vcf_1 <-  rTASSEL::readGenotypeTableFromPath(path = paste(vcfpath_goodman, "subset_goodman_snps_chr", 1, ".vcf.gz", sep = ""),
                                           keepDepth = FALSE)
print(vcf_1)

# # chrom 5, tissue LMAD
# vcf_5 <-  rTASSEL::readGenotypeTableFromPath(path = paste(vcfpath_goodman, "subset_goodman_snps_chr", 5, ".vcf.gz", sep = ""),
#                                              keepDepth = FALSE)
# print(vcf_5)


# # ------------
# # L3Base CHROM 1
# # ------------
# 
# # Join genotypes with (phenotypes + PCs + PEERs)
# tasPhenoDFGenotype_L3Base <- rTASSEL::readGenotypePhenotype(
#   genoPathOrObj = vcf_1,
#   phenoPathDFOrObj = tasPhenoDF_L3Base)
# print(tasPhenoDFGenotype_L3Base)
# 
# # Do a light MAF filter to remove invariant sites
# tasGenoPhenoFilt_L3Base <- rTASSEL::filterGenotypeTableSites(
#   tasObj = tasPhenoDFGenotype_L3Base,
#   siteRangeFilterType = "none",
#   siteMinAlleleFreq = 0.01,
#   siteMaxAlleleFreq = 1.0,
#   siteMinCount = 100)
# 
# # Run fast association, write files to disk.
# rTASSEL::assocModelFitter(
#   tasObj = tasGenoPhenoFilt_L3Base,
#   formula = . ~ ., 
#   fitMarkers = TRUE,
#   kinship = NULL,
#   fastAssociation = TRUE,
#   maxP = pValThresh,
#   maxThreads = nThreads,
#   outputFile = paste("chrom_", 1, "_fast_assoc_results_randomSNPs_", "L3Base", "_Kremling_2018", sep = ""))
# 
# 
# ---------------
# GShoot CHROM 1
# ---------------

# Join genotypes with (phenotypes + PCs + PEERs)
tasPhenoDFGenotype_GShoot <- rTASSEL::readGenotypePhenotype(
  genoPathOrObj = vcf_1,
  phenoPathDFOrObj = tasPhenoDF_GShoot)
print(tasPhenoDFGenotype_GShoot)

# Do a light MAF filter to remove invariant sites
tasGenoPhenoFilt_GShoot <- rTASSEL::filterGenotypeTableSites(
  tasObj = tasPhenoDFGenotype_GShoot,
  siteRangeFilterType = "none",
  siteMinAlleleFreq = 0.01,
  siteMaxAlleleFreq = 1.0,
  siteMinCount = 100)

# Run fast association, write files to disk.
rTASSEL::assocModelFitter(
  tasObj = tasGenoPhenoFilt_GShoot,
  formula = . ~ ., 
  fitMarkers = TRUE,
  kinship = NULL,
  fastAssociation = TRUE,
  maxP = pValThresh,
  maxThreads = nThreads,
  outputFile = paste("chrom_", 1, "_fast_assoc_results_randomSNPs_", "GShoot", "_Kremling_2018", sep = ""))

# 
# # ---------------
# # LMAD CHROM 1
# # ---------------
# 
# # Join genotypes with (phenotypes + PCs + PEERs)
# tasPhenoDFGenotype_LMAD <- rTASSEL::readGenotypePhenotype(
#   genoPathOrObj = vcf_1,
#   phenoPathDFOrObj = tasPhenoDF_LMAD)
# print(tasPhenoDFGenotype_LMAD)
# 
# # Do a light MAF filter to remove invariant sites
# tasGenoPhenoFilt_LMAD <- rTASSEL::filterGenotypeTableSites(
#   tasObj = tasPhenoDFGenotype_LMAD,
#   siteRangeFilterType = "none",
#   siteMinAlleleFreq = 0.01,
#   siteMaxAlleleFreq = 1.0,
#   siteMinCount = 100)
# 
# # Run fast association, write files to disk.
# rTASSEL::assocModelFitter(
#   tasObj = tasGenoPhenoFilt_LMAD,
#   formula = . ~ ., 
#   fitMarkers = TRUE,
#   kinship = NULL,
#   fastAssociation = TRUE,
#   maxP = pValThresh,
#   maxThreads = nThreads,
#   outputFile = paste("chrom_", 1, "_fast_assoc_results_randomSNPs_", "LMAD", "_Kremling_2018", sep = ""))
# 
# 
# # --------------
# # LMAD CHROM 5
# # --------------
# 
# # Join genotypes with (phenotypes + PCs + PEERs)
# tasPhenoDFGenotype_LMAD <- rTASSEL::readGenotypePhenotype(
#   genoPathOrObj = vcf_5,
#   phenoPathDFOrObj = tasPhenoDF_LMAD)
# print(tasPhenoDFGenotype_LMAD)
# 
# # Do a light MAF filter to remove invariant sites
# tasGenoPhenoFilt_LMAD <- rTASSEL::filterGenotypeTableSites(
#   tasObj = tasPhenoDFGenotype_LMAD,
#   siteRangeFilterType = "none",
#   siteMinAlleleFreq = 0.01,
#   siteMaxAlleleFreq = 1.0,
#   siteMinCount = 100)
# 
# # Run fast association, write files to disk.
# rTASSEL::assocModelFitter(
#   tasObj = tasGenoPhenoFilt_LMAD,
#   formula = . ~ ., 
#   fitMarkers = TRUE,
#   kinship = NULL,
#   fastAssociation = TRUE,
#   maxP = pValThresh,
#   maxThreads = nThreads,
#   outputFile = paste("chrom_", 5, "_fast_assoc_results_randomSNPs_", "LMAD", "_Kremling_2018", sep = ""))
# 
# 

# ----------------------------------------
# Process Results (to csv and rearranged)
# ----------------------------------------

# List all files in directory
expression_results <- "/workdir/mbb262/results/kremling/"
my_files <- list.files(expression_results, pattern = "[0-9]{4}.txt$")
my_files <- my_files[c(2)]

# Remove .txt extension
expression_short_names <- gsub("\\.txt$", "", my_files)

# expression result output directory
expression_output <- "/workdir/mbb262/results/kremling/processed_results/"

# Lapply loop to format results, output them to new directory
setwd("/workdir/mbb262/results/kremling/processed_results/")

lapply(seq_along(my_files),  function(i) {
  message("On trait: ", my_files[i])
  
  # Read in file
  expression_dat <- data.table::fread(paste0(expression_results, my_files[i]), nThread = 85) 
  
  # Make an end column
  expression_dat$end <- expression_dat$Pos
  
  # Rearrage columns and discard extras
  expression_dat <- expression_dat %>% select("Trait", "Chr", "Pos", "end")
  
  # Rename columns to work with Zack's pipeline
  colnames(expression_dat) <- c("trait", "seqid", "start", "end")
  
  # Export files to new directory
  data.table::fwrite(expression_dat, file = paste0(expression_output, expression_short_names[i], "_randomSNPs_reformatted_zack.csv"), nThread = 85)
})

# gzip then send back to blfs1
# gzip /workdir/mbb262/results/kremling/processed_results/*.csv
# scp /workdir/mbb262/results/kremling/processed_results/*.gz mbb262@cbfs1.biohpc.cornell.edu:/data1/users/mbb262/results/results_goodman_panel_kremling_maffilter/kremling_randomSNPs


# try to find problematic rows
# 143326Zm00001d049540_GShoot
lala1 <- temp[grep("143326Zm00001d049540_GShoot", temp$trait), ]
lala2 <- temp[grep("143326Zm00001d049540_GShoot", temp$seqid), ]
lala3 <- temp[grep("143326Zm00001d049540_GShoot", temp$start), ]
lala4 <- temp[grep("143326Zm00001d049540_GShoot", temp$end), ]


lala<- temp[grep("Zm00001d049540_GShoot", temp$trait), ]




# mapping two tissues, chrom 4 L3tip and kernel chrom 8
# something went wrong during data export and now won't load into R



# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2021-10-11
#
# Description 
#   - Preparing Kremiling data for fast association
#   - Running Kremling data in rTASSEL with the model:
#   - y ~ SNP + 3 gPCs (Ames to 282 method) + 25 PEER factors from Kremling et al.
#  Rscript /home/mbb262/git_projects/haplotype_v_snp_gwas/src/R/assoc_mapping/goodman282_gwas_rtassel_kremling.R > /workdir/mbb262/results_kremling/rtassel_kremling_goodman_r_output.txt
# ---------------------------------------------------------------

# Set directory
setwd("/workdir/mbb262/results/kremling_not_random/")

# Set memory limits
options(java.parameters = c("-Xmx500g"))

# Load in packages
library(DESeq2)
library(dplyr)
library(magrittr)
library(data.table)
library(rTASSEL)

# Start logging
rTASSEL::startLogger(fullPath = NULL, fileName = NULL)


# ----------------------------------
#   Load in pre-analyzed Covariates
# ----------------------------------

# Load in global PCs (Ames to Goodman method)
ames2goodman_gPCs_3gPCs <- read.csv("/home/mbb262/git_projects/pleiotropy/data/Q_all_3gPCs_allGoodman.csv", header = TRUE)
ames2goodman_gPCs_3gPCs <- ames2goodman_gPCs_3gPCs[,-c(5:6)]
colnames(ames2goodman_gPCs_3gPCs) <- c("Taxa", "PC1", "PC2", "PC3")

# Download in PEERs from irods
# icd /iplant/home/shared/commons_repo/curated/Kremling_Nature3RNASeq282_March2018/Expression_matrix/TASSEL_fmt_expression_w_covariates
# iget -r ./
# Remove first two lines with <Phenotype> from tassel header (I would load these into R with Rtassel but I keep getting an error)
# sed -i '1d' ./*.txt # ran twice because I'm not very clever
# sed -i '1d' ./*.txt

# Or, get PEERS from my blfs1 account --> do not have to remove first lines of file
# mkdir peers_copy
# scp -r mbb262@cbsublfs1.tc.cornell.edu:/data1/users/mbb262/phenotypes/kremling_2018_naturegenetics_expression/peers/*.txt /workdir/mbb262/peers_copy

# Load in peers
setwd("/workdir/mbb262/peers_copy/")
temp = list.files(pattern="*.txt")
myfiles = lapply(temp, 
                 data.table::fread, 
                 header=TRUE, 
                 sep = "\t", 
                 select = c(1, 7:31))

# Create shorter name for lists and then rename myfiles
short_names <- gsub("TASSEL_HEADER_df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames_", "", temp)# Merge the files together
short_names <- gsub("_sm_rand_and_box_coxed_w_MDSPCs_and_25PEERs_avgDups.txt", "", short_names)
names(myfiles) <- short_names

# Join with my PCs
myfiles <- lapply(myfiles, inner_join, x = ames2goodman_gPCs_3gPCs)


# ------------------------
# Format Kremling data
# ------------------------

# Code used to make these files is hashed out at the bottom of this script,
# these are the cleaned up files on blfsq
# scp -r mbb262@cbsublfs1.tc.cornell.edu:/data1/users/mbb262/phenotypes/kremling_2018_naturegenetics_expression/normalized_counts/kremling_raw_count_v4_hapmap321taxaid /workdir/mbb262

# Get all file names
setwd("/workdir/mbb262/phenotypes/kremling/")
temp = list.files("/workdir/mbb262/phenotypes/kremling/", pattern="*.csv")

# Load data into r as sepearte dfs within a single list
formatted_kremling <- lapply(temp, 
                             data.table::fread, 
                             header=TRUE)

# Shorten names
short_names <- gsub("_kremling_formatted_v4_hapmapids.csv", "", temp)
names(formatted_kremling) <- short_names # rename each element in the list to the tissue

# Remove L3 from the list
formatted_kremling <- formatted_kremling[-4]

# Add on tissue name to gene names
for (i in seq_along(formatted_kremling)){
  colnames(formatted_kremling[[i]]) <- paste0(colnames(formatted_kremling[[i]]), "_", names(formatted_kremling[i]))
  colnames(formatted_kremling[[i]])[1] <- "Taxa"
}

# Check
lapply(formatted_kremling, function(x) x[1:5,1:5])


# -------------------------------------
# Join expression data with PCs + PEERs
# -------------------------------------

# Check order of both lists, must be the same
names(formatted_kremling)
names(myfiles)

# Innerjoin two lists by index and by taxa name
expression_pcs_peers <- purrr::map2(myfiles, formatted_kremling, inner_join, by = "Taxa")

# Check
# str(expression_pcs_peers)
# names(expression_pcs_peers)
# expression_pcs_peers$GRoot[1:5,1:30]
# expression_pcs_peers$LMAN[1:5,1:30]
# str(expression_pcs_peers$LMAN[1:5,1:30])

# Tasselize merged phenotypes + PCs + PEERS
tasPhenoDF <- lapply(expression_pcs_peers, function(x) {rTASSEL::readPhenotypeFromDataFrame(
  phenotypeDF = x,
  taxaID = "Taxa",
  attributeTypes = c(rep("covariate", 28), rep("data", ncol(x)-29)))
})
tasPhenoDF
names(tasPhenoDF[1])

# Subset list to re-run problematic files (L3Tip, L3Base, GShoot, Kern)
tasPhenoDF <- list(tasPhenoDF$Kern, tasPhenoDF$L3Tip)
names(tasPhenoDF) <- c("Kern", "L3Tip")


# ---------------------------
#       Map phenotypes
# ---------------------------

# Local variables
# scp -r mbb262@cbsublfs1.tc.cornell.edu:/data1/users/mbb262/genotypes/goodman282/*_imputed_goodman282.vcf.gz /workdir/mbb262/genotypes/
vcf <- "/workdir/mbb262/genotypes/goodman/"

# Set directory for results
setwd("/workdir/mbb262/results/kremling_not_random/")

# Run fast association L3 Tip chrom 1
# Load in genotype table
vcf <-  rTASSEL::readGenotypeTableFromPath(path = paste(vcf,"hmp321_282_agpv4_merged_chr", 1, "_imputed_goodman282.vcf.gz", sep = ""),
                                           keepDepth = FALSE)
print(vcf)

# Iterate through phenotypes
j <- 2
message("I am on phenotype file ", names(tasPhenoDF)[j])

# Join genotypes with (phenotypes + PCs + PEERs)
tasPhenoDFGenotype <- rTASSEL::readGenotypePhenotype(
  genoPathOrObj = vcf,
  phenoPathDFOrObj = tasPhenoDF[[j]])
print(tasPhenoDFGenotype)

# Do a light MAF filter to remove invariant sites
tasGenoPhenoFilt <- rTASSEL::filterGenotypeTableSites(
  tasObj = tasPhenoDFGenotype,
  siteRangeFilterType = "none",
  siteMinAlleleFreq = 0.05,
  siteMaxAlleleFreq = 1.0,
  siteMinCount = 100)

# Run fast association, write files to disk.
rTASSEL::assocModelFitter(
  tasObj = tasGenoPhenoFilt,
  formula = . ~ ., 
  fitMarkers = TRUE,
  kinship = NULL,
  fastAssociation = TRUE,
  maxP = 0.001,
  maxThreads = 80,
  outputFile = paste("chrom_", 1, "_fast_assoc_results_", names(tasPhenoDF)[j], "_Kremling_2018", sep = ""))


# run fast association Kern chrom 8
# Load in genotype table
vcf <- "/workdir/mbb262/genotypes/goodman/"

# Set directory for results
setwd("/workdir/mbb262/results/kremling_not_random/")

vcf <-  rTASSEL::readGenotypeTableFromPath(path = paste(vcf,"hmp321_282_agpv4_merged_chr", 8, "_imputed_goodman282.vcf.gz", sep = ""),
                                           keepDepth = FALSE)
print(vcf)

j <- 1
message("I am on phenotype file ", names(tasPhenoDF)[j])

# Join genotypes with (phenotypes + PCs + PEERs)
tasPhenoDFGenotype <- rTASSEL::readGenotypePhenotype(
  genoPathOrObj = vcf,
  phenoPathDFOrObj = tasPhenoDF[[j]])
print(tasPhenoDFGenotype)

# Do a light MAF filter to remove invariant sites
tasGenoPhenoFilt <- rTASSEL::filterGenotypeTableSites(
  tasObj = tasPhenoDFGenotype,
  siteRangeFilterType = "none",
  siteMinAlleleFreq = 0.05,
  siteMaxAlleleFreq = 1.0,
  siteMinCount = 100)

# Run fast association, write files to disk.
rTASSEL::assocModelFitter(
  tasObj = tasGenoPhenoFilt,
  formula = . ~ ., 
  fitMarkers = TRUE,
  kinship = NULL,
  fastAssociation = TRUE,
  maxP = 0.001,
  maxThreads = 80,
  outputFile = paste("chrom_", 8, "_fast_assoc_results_", names(tasPhenoDF)[j], "_Kremling_2018", sep = ""))

# ----------------------------------------
# Process Results (to csv and rearranged)
# ----------------------------------------

# get from blfs1
# scp mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/mbb262/results/pleiotropy/gwa_results_goodman_panel_kremling_maffilter/unfiltered_goodman_kremling_results

# List all files in directory
expression_results <- "/workdir/mbb262/results/kremling_not_random/"
my_files <- list.files(expression_results, pattern = ".txt$")
my_files <- my_files[c(4,34)] # check to make sure I'm pulling out two correct files

# Remove .txt extension
expression_short_names <- gsub("\\.txt$", "", my_files)

# expression result output directory
expression_output <- "/workdir/mbb262/results/kremling_not_random/processed_results/"

# Lapply loop to format results, output them to new directory
setwd("/workdir/mbb262/results/kremling/processed_results/")

lapply(seq_along(my_files),  function(i) {
  message("On trait: ", my_files[i])
  
  # Read in file
  expression_dat <- data.table::fread(paste0(expression_results, my_files[i]), nThread = 85) 
  
  # Make an end column
  expression_dat$end <- expression_dat$Pos
  
  # Rearrage columns and discard extras
  expression_dat <- expression_dat %>% select("Trait", "Chr", "Pos", "end", "p")
  
  # Rename columns to work with Zack's pipeline
  colnames(expression_dat) <- c("trait", "seqid", "start", "end", "p")
  
  # Export files to new directory
  data.table::fwrite(expression_dat, file = paste0(expression_output, expression_short_names[i], "_reformatted_zack.csv"), nThread = 85)
})

