# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2021-06-09 
#
# Description 
#   - Example script for running fast association in rTassel
# ---------------------------------------------------------------

# Install package
if (!require("devtools")) install.packages("devtools")
if (!require("devtools")) install.packages("devtools")
devtools::install_bitbucket(
  repo = "bucklerlab/rtassel",
  build_vignettes = F
)

# Setting memory BEFORE loading in rtassel
options(java.parameters = c("-Xmx250g"))

# Load package
library(rTASSEL)
library(dplyr)

# Start logging
rTASSEL::startLogger(fullPath = "/workdir/mbb262/results_goodman_panel/", 
                     fileName = NULL)


# ----------------------------------
#   Load in pre-analyzed Covariates
# ----------------------------------

# Load in pre-analyzed PCs
# 3-5 PCs in the 282 works great


# ---------------------
# Analysis Pipeline
# ---------------------

# For-loop iterates over chromosomes
# Loads in phenotype files
# Intersects phenotypes, genotypes, global PCs, window PCs
# Runs through fast association and saves results to file

# Vcf path
vcfpath <- "/workdir/mbb262/goodman282/"
phenopath <- "/workdir/mbb262/tasselized_for_fa/"

# List of metabolites to analyze 
# (each of these files contain the 282 genotypes by 700-1000 different metabolites, you can just have 1 file
# with many traits if that works for you) 
meta <- c("15_goodman_phenotypes_raw_posLeafTip_Zhou_2019.txt", 
          "16_goodman_phenotypes_negLeafBase_Zhou_2019.txt", 
          "17_goodman_phenotypes_posLeafBase_Zhou_2019.txt",
          "14_goodman_phenotypes_raw_negLeafTip_Zhou_2019.txt")


# Set directory to where results will be output
setwd("/workdir/mbb262/results_goodman_panel") # Where to output results

# For loop to iterate through all chromosomes then all phenotypes to run fast association
# Analyze metabolites 
for (i in seq_len(10)){
  message("I am on chromosome ", i)
  
  # Load in genotype table
  vcf <-  rTASSEL::readGenotypeTableFromPath(path = paste(vcfpath,"hmp321_282_agpv4_merged_chr", i, "_imputed_goodman282.vcf.gz", sep = ""),
                                             keepDepth = FALSE)
  print(vcf)
  
  # Iterate through phenotypes
  for (j in seq_along(meta)){
    message("I am on phenotype file ", meta[j])
    
    # Load in phenotype file
    tas_phenotype <- rTASSEL::readPhenotypeFromPath(path = paste(phenopath, meta[j], sep = ""))
    print(tas_phenotype)
    
    # Turn into r dataframe (code needs to be updated, use hacky way for now)
    # phenotype <- rTASSEL::getPhenotypeDF(tasObj = tas_phenotype)
    phenotype <- data.frame(rTASSEL:::tableReportToDF(tas_phenotype@jPhenotypeTable))
    
    # Intersect phenotypes and PCs
    pheno_gpcs_windowPCs <- merge(x = phenotype, y = ames2goodman_gPCs_3gPCs)
    print(str(phenotype))
    
    # Tasselize merged phenotypes + global PCs + window PCs
    # attributeTypes assigns what is data and what is a covariate
    tasPhenoDF <- rTASSEL::readPhenotypeFromDataFrame(
      phenotypeDF = pheno_gpcs_windowPCs,
      taxaID = "Taxa",
      attributeTypes = c(rep("data", ncol(phenotype)-1), rep("covariate", ncol(ames2goodman_gPCs_3gPCs)-1)))
    
    # Join genotypes with (phenotypes + g PCs)
    tasPhenoDF <- rTASSEL::readGenotypePhenotype(
      genoPathOrObj = vcf,
      phenoPathDFOrObj = tasPhenoDF)
    print(tasPhenoDF)
    
    # Do a light MAF filter to remove invariant sites
    tasGenoPhenoFilt <- rTASSEL::filterGenotypeTableSites(
      tasObj = tasPhenoDF,
      siteRangeFilterType = "none",
      siteMinAlleleFreq = 0.01,
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
      maxThreads = 39,
      outputFile = paste("chrom_", i, "_fast_assoc_results_", meta[j], sep = "")
    )
  }
}
