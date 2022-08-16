# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-08-15
# Updated... 2022-08-15

# Description:
# Run a mixed model in rTASSEL
# phenotype ~ SNP + Kinship
# ---------------------------------------------------------------


# Setting memory
options(java.parameters = c("-Xmx120g"))
numThreads <- 5
p_value <- 0.01

# Load package
library(rTASSEL)
library(dplyr)
library(ggplot2)

# Start logging
rTASSEL::startLogger(fullPath = NULL, fileName = NULL)


# ------------------------------------------------
# Gather phenotypes, genotypes
# ------------------------------------------------

# Load in a numeric genotype table
numericpath <- "~/Desktop/TutorialData/numeric_geno.txt"
genotype <-  rTASSEL::readGenotypeTableFromPath(path = numericpath, keepDepth = FALSE)

# Load in phenotype file
phenotype <- data.table::fread("/workdir/mbb262/Mann_GWA/Pheno-Absolute-log10.csv", nThread = numThreads, drop = c(1))
colnames(phenotype) <- c("Taxa", "pheno")

# Tasselize phenotype
tasPhenoDF <- rTASSEL::readPhenotypeFromDataFrame(
  phenotypeDF = phenotype,
  taxaID = "Taxa",
  attributeTypes = c("data"))


# ----------------------------------------------------------------------------------------------
# Run Mixed linear model
# phenotype ~ SNP + Kinship
# ----------------------------------------------------------------------------------------------

# Join genotypes with phenotypes
tasPhenoDF <- rTASSEL::readGenotypePhenotype(
  genoPathOrObj = genotype,
  phenoPathDFOrObj = tasPhenoDF)

# Do a light MAF filter to remove invariant sites
# tasGenoPhenoFilt <- rTASSEL::filterGenotypeTableSites(
#   tasObj = tasPhenoDF,
#   siteRangeFilterType = "none",
#   siteMinAlleleFreq = 0.01,
#   siteMaxAlleleFreq = 1.0,
#   siteMinCount = 13)
print(tasGenoPhenoFilt)

# Calculate kinship matrix
tasKin <- kinshipMatrix(tasObj = tasGenoPhenoFilt)

# Run fast association, write files to disk.
rTASSEL::assocModelFitter(
  tasObj = tasGenoPhenoFilt,
  formula = pheno ~ .,
  fitMarkers = TRUE,
  kinship = tasKin,
  fastAssociation = TRUE,
  maxP = p_value,
  maxThreads = numThreads,
  outputFile = "/workdir/mbb262/fast_assoc_results_model_5.txt")

# Load in association results
gwas_result_fa_model_5 <- data.table::fread("/workdir/mbb262/fast_assoc_results_model_5.txt")
