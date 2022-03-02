phenos_nam <- c("01_nam_phenotypes_Buckler_2009.txt")

# Path to files
phenopath_nam <- "~/Downloads/workdir/mbb262/phenotypes/nam/"

# Path to PCs
global_pcs <- data.table::fread("~/Downloads/workdir/mbb262/03_pcs/nam/ames2nam_3gPCs_allNAM_gPCs_only_not4tassel.txt")
main_pcs <- data.table::fread("~/Downloads/workdir/mbb262/03_pcs/nam/mainGeneWindow360_ames2nam_local_3allNAM_varaibleNumPCs_not4tassel.txt")
midway_pcs <- data.table::fread("~/Downloads/workdir/mbb262/03_pcs/nam/midwayGeneWindow360_ames2nam_local_3allNAM_varaibleNumPCs_not4tassel.txt")

# variables
j <- 1
p <- 1

# Load in phenotype file 
phenotype <- data.table::fread(paste(phenopath_nam, phenos_nam[j], sep = ""))
phenotype <- phenotype[,1:2]

# Change taxa name (currently is <Trait>)
colnames(phenotype)[1] <- "Taxa"

# join phenotypes and PCs
phenotype <- merge(phenotype, global_pcs, by = "Taxa")
pheno_main_pcs <- merge(phenotype, main_pcs, by = "Taxa")
pheno_midway_pcs <- merge(phenotype, midway_pcs, by = "Taxa")

# Fit two models to the same trait??
fit_main <- lm(DTA_BLUP.Buckler_2009_nam ~ ., data = pheno_main_pcs[,-1])
fit_midway <- lm(DTA_BLUP.Buckler_2009_nam ~ ., data = pheno_midway_pcs[,-1])

# Grab fitted values and residuals for both models??
y_hat_main <- cbind(pheno_main_pcs[,1], fit_main$fitted.values)
residual_y_main <- cbind(pheno_main_pcs[,1], fit_main$residuals)
y_hat_midway <- cbind(pheno_midway_pcs[,1], fit_midway$fitted.values)
residual_y_midway <- cbind(pheno_midway_pcs[,1], fit_midway$residuals)

# Change column names to traits
colnames(y_hat)[2] <- trait_ids[p]
colnames(residual_y)[2] <- trait_ids[p]


# Permute trait 10X times
temp <- lapply(seq_len(10), function(perm) {
  message("Permutation: ", perm)
  
  # Collect phenotypes outside of loop
  perm_adjusted_y_hat <- pcs_pheno[,1] %>% data.frame()
  
  # subsample rows
  perm_df <- cbind(residual_y[,1], residual_y[sample(nrow(residual_y)), -1])
  
  # Add back y_hat
  add_y_hat <- merge(perm_df, y_hat, by = "Taxa")
  add_y_hat <- cbind(add_y_hat[,1], rowSums(add_y_hat[,-1]))
  
  # Change column names
  colnames(add_y_hat)[2] <- paste0(trait_ids[p], "_permutation_", perm)
  
  # merge into outside list
  perm_adjusted_y_hat <- merge(perm_adjusted_y_hat, add_y_hat, by = "taxa")
  return(perm_adjusted_y_hat)
})

# Combine all permutations for a single trait
lala <- Reduce(function(x, y) merge(x, y, by = "taxa") , temp)