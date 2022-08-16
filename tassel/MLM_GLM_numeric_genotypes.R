# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-08-15
# Updated... 2022-08-15

# Description:
# Run a mixed model in rTASSEL
# phenotype ~ numeric genotype PAV + Kinship
# 
# Run GLM
# phenotype ~ numeric genotype PAV
#
# rtassel wiki
# https://maize-genetics.github.io/rTASSEL/articles/rtassel_walkthrough.html
# ---------------------------------------------------------------


# Setting memory
options(java.parameters = c("-Xmx30g"))
numThreads <- 5
p_value <- 0.01

# Load package
library(rTASSEL)

# Start logging
rTASSEL::startLogger(fullPath = NULL, fileName = NULL)


# ------------------------------------------------
# Gather phenotypes, genotypes
# ------------------------------------------------

# Load in a numeric genotype table
genotype <-  rTASSEL::readGenotypeTableFromPath(path = "/workdir/mbb262/mdp_genotype_with_Probability.txt")

# Load in phenotype file
phenotype <- data.table::fread("/workdir/mbb262/blued_phenotypes.txt")
phenotype <- phenotype[,c(1:2)] # kust keep taxa names and 1 phenotype

# Tasselize phenotype
tasPhenoDF <- rTASSEL::readPhenotypeFromDataFrame(
  phenotypeDF = phenotype,
  taxaID = "Taxa",
  attributeTypes = c("data"))

# Read in kinship matrix
tasKin <- readTasselDistanceMatrix("/workdir/mbb262/mdp_kinship.txt")


# ----------------------------------------------------------------------------------------------
# Run Mixed linear model
# phenotype ~ SNP + Kinship
# ----------------------------------------------------------------------------------------------

# Join genotypes with phenotypes
tasPhenoDF <- rTASSEL::readGenotypePhenotype(
  genoPathOrObj = genotype,
  phenoPathDFOrObj = tasPhenoDF)
print(tasPhenoDF)

# Run MLM, write files to disk.
rTASSEL::assocModelFitter(
  tasObj = tasPhenoDF,
  formula = EarDia ~ .,
  fitMarkers = TRUE,
  kinship = tasKin,
  fastAssociation = FALSE,
  maxP = p_value,
  maxThreads = numThreads,
  outputFile = "/workdir/mbb262/model_results.txt")


# Load in association results
gwas_results <- data.table::fread("/workdir/mbb262/model_results.txt")


# Calculate GLM, save to r environment
tasGLM <- rTASSEL::assocModelFitter(
  tasObj = tasPhenoDF,
  formula = EarDia ~ .,
  fitMarkers = TRUE, 
  kinship = NULL,
  fastAssociation = FALSE
)
