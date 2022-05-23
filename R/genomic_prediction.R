# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-01-24 
#
# Description 
# Genomic prediction using rTASSEL
#
# Adapted from:
# https://maize-genetics.github.io/rTASSEL/articles/rtassel_walkthrough.html#genomic-prediction
# ---------------------------------------------------------------

## Load packages
library("rtassel")


# ----------------------------
## Load phenotypes
# ----------------------------



# ----------------------------
## Load in genotypes
# ----------------------------

# Load in genotypes
tasGenoPheno

# Filter genotypes on MAF
tasGenoPhenoFilt <- rTASSEL::filterGenotypeTableSites(
  tasObj = tasGenoPheno,
  siteMinAlleleFreq = 0.05,
  siteMaxAlleleFreq = 1.0,
  siteRangeFilterType = "none"
)
tasGenoPhenoFilt

# Filter data on unique bed file
tasGenoPheno_bed <- rTASSEL::filterGenotypeTableSites(
  tasObj = tasGenoPhenoFilt,
  siteMinAlleleFreq = 0.05,
  siteMaxAlleleFreq = 1.0,
  siteRangeFilterType = "none" # replace with bed file
)
tasGenoPheno_bed


# ----------------------------
## Create kinship matrix
# ----------------------------

# Genome wide
k_genome <- kinshipMatrix(tasObj = tasGenoPhenoFilt)

# Subsetted genotypes
k_genome_subset <- kinshipMatrix(tasObj = tasGenoPheno_bed)

# Correlate kinship matrices
cor_matrices <- cor(c(k_genome), c(k_genome_subset))


# ----------------------------
## Run genomic prediction
# ----------------------------

# Results show correlation (r) of the observed values and predicted values is 
# calculated for the validation set. The mean and standard deviation of the mean 
# of the r’s are calculated for each trait and reported in the comments section 
# of the “Accuracy” data set that is output by the analysis. 
tasCV <- genomicPrediction(
  tasPhenoObj = tasGenoPheno,
  kinship     = tasKin,
  doCV        = TRUE,
  kFolds      = 5,
  nIter       = 1
)
head(tasCV)


# ----------------------------
## Plot results
# ----------------------------

