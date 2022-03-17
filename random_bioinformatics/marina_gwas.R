# Setting memory
options(java.parameters = c("-Xmx120g"))
numThreads <- 20
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

# Load in genotype table
vcfpath <- "/workdir/mbb262/Mann_GWA/"
vcf <-  rTASSEL::readGenotypeTableFromPath(path = paste(vcfpath,"DC_BEAGLE_OUTPUT.vcf.gz", sep = ""), keepDepth = FALSE)

# Load in phenotype file
phenotype <- data.table::fread("/workdir/mbb262/Mann_GWA/Pheno-Absolute-log10.csv", nThread = numThreads, drop = c(1))
colnames(phenotype) <- c("Taxa", "pheno")

# Load in pcs
pc_file <- data.table::fread("/workdir/mbb262/Mann_GWA/pc_file.txt")
colnames(pc_file)[1] <- "Taxa"
pc_file <- pc_file[,1:6] # Only keep the first 5 PCs


# ------------------------------------------------
# Run Fast Association
# y ~ SNP + PC1 + PC2 
# ------------------------------------------------

# Intersect phenotypes and PCs
pheno_gpcs <- merge(x = phenotype, y = pc_file)

# Tasselize merged phenotypes + global PCs 
tasPhenoDF <- rTASSEL::readPhenotypeFromDataFrame(
  phenotypeDF = pheno_gpcs,
  taxaID = "Taxa",
  attributeTypes = c(rep("data", 1), rep("covariate", 5)))

# Join genotypes with (phenotypes + g PCs)
tasPhenoDF <- rTASSEL::readGenotypePhenotype(
  genoPathOrObj = vcf,
  phenoPathDFOrObj = tasPhenoDF)

# Do a light MAF filter to remove invariant sites
tasGenoPhenoFilt <- rTASSEL::filterGenotypeTableSites(
  tasObj = tasPhenoDF,
  siteRangeFilterType = "none",
  siteMinAlleleFreq = 0.01,
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


# ------------------------------------------------
# Run Mixed linear model
# y ~ SNP + PC1 + PC2 + management + SNP*management
# ------------------------------------------------

# Tasselize merged phenotypes + global PCs 
tasPhenoDF <- rTASSEL::readPhenotypeFromDataFrame(
  phenotypeDF = phenotype,
  taxaID = "Taxa",
  attributeTypes = c("data", "factor"))

# Join genotypes with (phenotypes + g PCs)
tasPhenoDF <- rTASSEL::readGenotypePhenotype(
  genoPathOrObj = vcf,
  phenoPathDFOrObj = tasPhenoDF)

# Do a light MAF filter to remove invariant sites
tasGenoPhenoFilt <- rTASSEL::filterGenotypeTableSites(
  tasObj = tasPhenoDF,
  siteRangeFilterType = "none",
  siteMinAlleleFreq = 0.01,
  siteMaxAlleleFreq = 1.0,
  siteMinCount = 13)
print(tasGenoPhenoFilt)

# Calculate kinship matrix
tasKin <- kinshipMatrix(tasObj = tasGenoPhenoFilt)

# Calculate MLM
tasMLM <- rTASSEL::assocModelFitter(
  tasObj = tasGenoPhenoFilt,
  formula = pheno ~ management + PC1 + PC2,
  fitMarkers = TRUE,                 # <- set this to TRUE for GLM
  kinship = tasKin,
  fastAssociation = FALSE
)

# Load in association resuilts
gwas_result <- data.table::fread("/workdir/mbb262/mlm_results.txt")


# ----------------------------------
# Prepare to plot things
# ----------------------------------

# Format x axis to plot by chromosome
data_cum <- gwas_result %>% 
  group_by(Chr) %>% 
  summarise(max_bp = max(Pos)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  select(Chr, bp_add)

gwas_result <- gwas_result %>% 
  inner_join(data_cum, by = "Chr") %>% 
  mutate(bp_cum = Pos + bp_add)
  
# Plot everything out
ggplot(gwas_result, aes(x = bp_cum, y = -log10(p), color = as.factor(Chr), size = -log10(p))) +
    geom_point(size = 1) +
    geom_hline(yintercept = 5, color = "grey40", linetype = "dashed") +
    labs(x = NULL, y = "-log10(p)", title = "Marina GWA") + 
    theme_minimal()
  
  
  
  

