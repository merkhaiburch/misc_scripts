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
pc_file <- pc_file[,1:3] # Only keep the first 5 PCs

# Management file
management <- data.table::fread("/workdir/mbb262/Mann_GWA/env.txt", drop = c(1))
colnames(management) <- c("Taxa", "field")


# ------------------------------------------------
# Run Fast Association
# Model 1: y ~ SNP + management
# ------------------------------------------------

# Intersect phenotypes and PCs
pheno_gpcs <- merge(x = phenotype, y = management)
pheno_gpcs$field <- as.numeric(as.factor(pheno_gpcs$field))

# Tasselize merged phenotypes + global PCs 
tasPhenoDF <- rTASSEL::readPhenotypeFromDataFrame(
  phenotypeDF = pheno_gpcs,
  taxaID = "Taxa",
  attributeTypes = c("data", "covariate"))

# Join genotypes with (phenotypes + g PCs)
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
  outputFile = "/workdir/mbb262/fast_assoc_results_model_1.txt")

# Load in association resuilts
gwas_result_fa_model_1 <- data.table::fread("/workdir/mbb262/fast_assoc_results_model_1.txt")

# Format x axis to plot by chromosome
data_cum <- gwas_result_fa_model_1 %>% 
  group_by(Chr) %>% 
  summarise(max_bp = max(Pos)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  select(Chr, bp_add)

gwas_result <- gwas_result_fa_model_1 %>% 
  inner_join(data_cum, by = "Chr") %>% 
  mutate(bp_cum = Pos + bp_add)

# Plot everything out
model1 <- ggplot(gwas_result, aes(x = bp_cum, y = -log10(p), color = as.factor(Chr), size = -log10(p))) +
  geom_point(size = 1) +
  geom_hline(yintercept = 5, color = "grey40", linetype = "dashed") +
  labs(x = NULL, y = "-log10(p)", title = "y ~ SNP + management") + 
  theme_minimal()


# ------------------------------------------------
# Run Fast Association
# Model 2: y ~ SNP + PC1 + PC2 
# ------------------------------------------------

# Intersect phenotypes and PCs
pheno_gpcs <- merge(x = phenotype, y = pc_file)

# Tasselize merged phenotypes + global PCs 
tasPhenoDF <- rTASSEL::readPhenotypeFromDataFrame(
  phenotypeDF = pheno_gpcs,
  taxaID = "Taxa",
  attributeTypes = c("data", "covariate","covariate"))

# Join genotypes with (phenotypes + g PCs)
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

# Load in association resuilts
gwas_result_fa_model_2 <- data.table::fread("/workdir/mbb262/fast_assoc_results_model_2.txt")

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


# ------------------------------------------------
# Run Fast Association
# Model 3: y ~ SNP + PC1 + PC2 + management
# ------------------------------------------------

# Intersect phenotypes and PCs
pheno_gpcs <- merge(x = phenotype, y = pc_file)
pheno_gpcs <- merge(x = pheno_gpcs, y = management)
pheno_gpcs$field <- as.numeric(as.factor(pheno_gpcs$field))

# Tasselize merged phenotypes + global PCs 
tasPhenoDF <- rTASSEL::readPhenotypeFromDataFrame(
  phenotypeDF = pheno_gpcs,
  taxaID = "Taxa",
  attributeTypes = c("data", "covariate","covariate", "covariate"))

# Join genotypes with (phenotypes + g PCs)
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
  outputFile = "/workdir/mbb262/fast_assoc_results_model_3.txt")

# Load in association resuilts
gwas_result_fa_model_3 <- data.table::fread("/workdir/mbb262/fast_assoc_results_model_3.txt")

# Format x axis to plot by chromosome
data_cum <- gwas_result_fa_model_3 %>% 
  group_by(Chr) %>% 
  summarise(max_bp = max(Pos)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  select(Chr, bp_add)

gwas_result <- gwas_result_fa_model_3 %>% 
  inner_join(data_cum, by = "Chr") %>% 
  mutate(bp_cum = Pos + bp_add)

# Plot everything out
model3 <- ggplot(gwas_result, aes(x = bp_cum, y = -log10(p), color = as.factor(Chr), size = -log10(p))) +
  geom_point(size = 1) +
  geom_hline(yintercept = 5, color = "grey40", linetype = "dashed") +
  labs(x = NULL, y = "-log10(p)", title = "y ~ SNP + PC1 + PC2 + management") + 
  theme_minimal()


# ------------------------------------------------
# Run Fast Association
# Model 4: BLUE ~ SNP + PC1 + PC2
# ------------------------------------------------

# Intersect phenotypes and PCs
pheno_manage <- merge(x = phenotype, y = management)

# Specify data types
pheno_manage <- rTASSEL::readPhenotypeFromDataFrame(
  phenotypeDF = pheno_manage,
  taxaID = "Taxa",
  attributeTypes = c("data", "factor"))

# Calculate BLUEs
tasBLUE <- rTASSEL::assocModelFitter(
  tasObj = pheno_manage,
  formula = . ~ .,                  # <- All data is used!
  fitMarkers = FALSE,
  kinship = NULL,
  fastAssociation = FALSE
)

# Merge blues with pcs
tasPhenoDF <- merge(x = tasBLUE$BLUE, y = pc_file)

tasPhenoDF <- rTASSEL::readPhenotypeFromDataFrame(
  phenotypeDF = tasPhenoDF,
  taxaID = "Taxa",
  attributeTypes = c("data", "covariate", "covariate"))

# Join genotypes with (phenotypes + g PCs)
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
  outputFile = "/workdir/mbb262/fast_assoc_results_model_4.txt")

# Load in association resuilts
gwas_result_fa_model_4 <- data.table::fread("/workdir/mbb262/fast_assoc_results_model_4.txt")

# Format x axis to plot by chromosome
data_cum <- gwas_result_fa_model_4 %>% 
  group_by(Chr) %>% 
  summarise(max_bp = max(Pos)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  select(Chr, bp_add)

gwas_result <- gwas_result_fa_model_4 %>% 
  inner_join(data_cum, by = "Chr") %>% 
  mutate(bp_cum = Pos + bp_add)

# Plot everything out
model4 <- ggplot(gwas_result, aes(x = bp_cum, y = -log10(p), color = as.factor(Chr), size = -log10(p))) +
  geom_point(size = 1) +
  geom_hline(yintercept = 5, color = "grey40", linetype = "dashed") +
  labs(x = NULL, y = "-log10(p)", title = "BLUE ~ SNP + PC1 + PC2") + 
  theme_minimal()


# Plot all things together
model1/model2/model3/model4+ 
  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom')

ggpubr::ggarrange(model1, model2, model3, model4, nrow = 4, ncol = 1,  common.legend = TRUE)



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