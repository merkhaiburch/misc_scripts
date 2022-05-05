# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2021-03-25
# Updated... 2022-03-25
#
# Description 
# Run genome wide association using rTassel
# rtassel website: https://github.com/maize-genetics/rTASSEL
# tutorial: https://maize-genetics.github.io/rTASSEL/articles/rtassel_walkthrough.html
# ---------------------------------------------------------------

# Setting memory
options(java.parameters = c("-Xmx120g"))
numThreads <- 20
p_value <- 0.00001

# Load package
library(rTASSEL)
library(dplyr)
library(ggplot2)

# Start logging
rTASSEL::startLogger(fullPath = NULL, fileName = NULL)


# ------------------------------------------------
# Gather phenotypes, genotypes
# ------------------------------------------------

# Load in genotype table
vcfpath <- "/path/to/dir/"
vcf <-  rTASSEL::readGenotypeTableFromPath(path = paste(vcfpath,"vcf.vcf.gz", sep = ""), 
                                           keepDepth = FALSE)

# Load in phenotype file
phenotype <- data.table::fread("/path/pheno.csv", nThread = numThreads)
colnames(phenotype) <- c("Taxa", "pheno")

# Load in pcs
pc_file <- data.table::fread("/path/pc_file.txt")
colnames(pc_file)[1] <- "Taxa"
pc_file <- pc_file[,1:3] # Only keep the first 3 PCs

# Intersect phenotypes and PCs
pheno_gpcs <- merge(x = phenotype, y = pc_file)


# ------------------------------------------------
# Run Fast Association
# Model: y ~ SNP + PC1 + PC2 
# ------------------------------------------------

# Tasselize merged phenotypes + global PCs 
tasPhenoDF <- rTASSEL::readPhenotypeFromDataFrame(
  phenotypeDF = pheno_gpcs,
  taxaID = "Taxa",
  attributeTypes = c("data", "covariate","covariate"))

# Join genotypes with (phenotypes + PCs)
tasPhenoDF <- rTASSEL::readGenotypePhenotype(
  genoPathOrObj = vcf,
  phenoPathDFOrObj = tasPhenoDF)

# Do a light MAF filter to remove invariant sites
tasGenoPhenoFilt <- rTASSEL::filterGenotypeTableSites(
  tasObj = tasPhenoDF,
  siteRangeFilterType = "none",
  siteMinAlleleFreq = 0.1,
  siteMaxAlleleFreq = 1.0,
  siteMinCount = 13)
print(tasGenoPhenoFilt)

# Run fast association, write files to disk.
rTASSEL::assocModelFitter(
  tasObj = tasGenoPhenoFilt,
  formula = pheno ~ .,
  fitMarkers = TRUE,
  kinship = NULL,
  fastAssociation = TRUE,
  maxP = p_value,
  maxThreads = numThreads,
  outputFile = "/workdir/mbb262/fast_assoc_results.txt")


# ----------------------------------------------
#           Plot Manhattan plot
# ----------------------------------------------

# Load in association resuilts
gwas_result_fa_model_2 <- data.table::fread("/workdir/mbb262/fast_assoc_results.txt")

# Format x axis to plot by chromosome
data_cum <- gwas_result_fa_model_2 %>% 
  group_by(Chr) %>% 
  summarise(max_bp = max(Pos)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  select(Chr, bp_add)

gwas_result <- gwas_result_fa_model_2 %>% 
  inner_join(data_cum, by = "Chr") %>% 
  mutate(bp_cum = Pos + bp_add)

# Plot everything out
model2 <- ggplot(gwas_result, aes(x = bp_cum, y = -log10(p), color = as.factor(Chr), size = -log10(p))) +
  geom_point(size = 1) +
  geom_hline(yintercept = 5, color = "grey40", linetype = "dashed") +
  labs(x = NULL, y = "-log10(p)", title = "y ~ SNP + PC1 + PC2") + 
  theme_minimal()
