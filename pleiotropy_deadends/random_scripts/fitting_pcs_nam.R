phenos_nam <- c("01_nam_phenotypes_Buckler_2009.txt", 
                "02_nam_phenotypes_Peiffer_2013.txt",
                "03_nam_phenotypes_Hung_2012_1.txt", 
                "04_nam_phenotypes_Wallace_2014.txt",
                "05_nam_phenotypes_Peiffer_2014.txt", 
                "06_nam_phenotypes_Hung_2012.txt",
                "07_nam_phenotypes_Poland_2011.txt", 
                "08_nam_phenotypes_Brown_2011.txt",
                "09_nam_phenotypes_Bian_2014_Kump_2011.txt", 
                "10_nam_phenotypes_Cook_2012.txt",
                "11_nam_phenotypes_Olukolu_2014.txt", 
                "12_nam_phenotypes_Benson_2015.txt",
                "13_nam_phenotypes_Diepenbrock_2017.txt", 
                "14_nam_phenotypes_Tian_2011.txt", 
                "15_nam_phenotypes_Olukolu_2016.txt", 
                "16_nam_phenotypes_Foerster_2015.txt")

# Format names for export
phenos_nam_export <- gsub("phenotypes_", "phenotypes_permuted10x_", phenos_nam)

# Path to files
phenopath_nam <- "~/Downloads/workdir/mbb262/phenotypes/nam/"

# Path to PCs
global_pcs <- data.table::fread("~/Downloads/workdir/mbb262/03_pcs/nam/ames2nam_3gPCs_allNAM_gPCs_only_not4tassel.txt")
main_pcs <- data.table::fread("~/Downloads/workdir/mbb262/03_pcs/nam/mainGeneWindow360_ames2nam_local_3allNAM_varaibleNumPCs_not4tassel.txt")
midway_pcs <- data.table::fread("~/Downloads/workdir/mbb262/03_pcs/nam/midwayGeneWindow360_ames2nam_local_3allNAM_varaibleNumPCs_not4tassel.txt")

# variables
j <- 1
p <- 1


# Iterate through phenotypes
pheno_files <- lapply(seq_len(length(phenos_nam)), function(j) {
  message("I am on phenotype file ", phenos_nam[j])
  
  # Load in phenotype file 
  phenotype <- data.table::fread(paste(phenopath_nam, phenos_nam[j], sep = ""))
  phenotype <- phenotype[,1:3] # SUBSET JUST FOR TESTING
  
  # Swap out semicolons in name with underscors
  colnames(phenotype) <- gsub(";", "_", colnames(phenotype))
  
  # Change taxa name (currently is <Trait>)
  colnames(phenotype)[1] <- "Taxa"
  
  # join phenotypes and PCs
  phenotype_gpcs <- merge(phenotype, global_pcs, by = "Taxa")
  pheno_main_pcs <- merge(phenotype_gpcs, main_pcs, by = "Taxa")
  pheno_midway_pcs <- merge(phenotype_gpcs, midway_pcs, by = "Taxa")
  
  # Gather trait IDs
  trait_ids <- colnames(pheno_main_pcs[, 2:ncol(phenotype)])
  
  # Iterate through each phenotype file, then each trait, permute 10x
  all_traits_perms <- lapply(seq_len(length(trait_ids)), function(p) {
    message("Permuting and adjusting trait: ", trait_ids[p])
    fit_main <- lm(paste(trait_ids[p], "~ ."), data = pheno_main_pcs[,-1]) #-1 removes Taxa ids
    fit_midway <- lm(paste(trait_ids[p], "~ ."), data = pheno_midway_pcs[,-1]) #-1 removes Taxa ids
    
    # Grab fitted values and residual
    y_hat_main <- cbind(pheno_main_pcs[,1], fit_main$fitted.values)
    residual_y_main <- cbind(pheno_main_pcs[,1], fit_midway$residuals)
    y_hat_midway <- cbind(pheno_midway_pcs[,1], fit_midway$fitted.values)
    residual_y_midway <- cbind(pheno_midway_pcs[,1], fit_midway$residuals)
    
    # Change column names to traits
    colnames(y_hat_main)[2] <- trait_ids[p]
    colnames(residual_y_main)[2] <- trait_ids[p]
    colnames(y_hat_midway)[2] <- trait_ids[p]
    colnames(residual_y_midway)[2] <- trait_ids[p]
    
    # Permute trait 10X times
    temp <- lapply(seq_len(10), function(perm) {
      message("Permutation: ", perm)
      
      # Create an order to permute traits with
      sample_order <- data.frame("Taxa" = sample(residual_y_main$Taxa))
      
      # Collect phenotypes
      perm_adjusted_y_hat_main <- pheno_main_pcs[,1] %>% data.frame()
      perm_adjusted_y_hat_midway <- pheno_midway_pcs[,1] %>% data.frame()
      
      # subsample rows in the same order between main and midway
      residual_y_main <- residual_y_main[order(match(residual_y_main$Taxa, sample_order$Taxa)), ]
      residual_y_midway <- residual_y_midway[order(match(residual_y_midway$Taxa, sample_order$Taxa)), ]
      
      # Add back y_hat
      add_y_hat_main <- merge(perm_df, y_hat, by = "Taxa")
      add_y_hat_main <- cbind(add_y_hat_main[,1], rowSums(add_y_hat_main[,-1]))
      
      add_y_hat <- merge(perm_df, y_hat, by = "Taxa")
      add_y_hat <- cbind(add_y_hat[,1], rowSums(add_y_hat[,-1]))
      
      # Change column names
      colnames(add_y_hat)[2] <- paste0(trait_ids[p], "_permutation_", perm)
      
      # merge into outside list
      perm_adjusted_y_hat <- merge(perm_adjusted_y_hat, add_y_hat, by = "taxa")
      # return(perm_adjusted_y_hat)
    })
    # Combine all permutations for a single trait
    lala <- Reduce(function(x, y) merge(x, y, by = "taxa") , temp)
  })
  # Combine all traits and their corresponding permutations
  temp2 <- Reduce(function(x, y) merge(x, y, by = "taxa") , all_traits_perms)
  
  # Export file
  data.table::fwrite(temp2,
                     paste0("~/Downloads/workdir/mbb262/phenotypes/02_permuted_fast_association_data/goodman/", phenos_good_export[j]))
  
})
  


# ----------------------------------------------------------
# Proof of concept reordering used in above function works
# ----------------------------------------------------------

test1 <- data.frame("Taxa" = c("A", "B", "C", "D"), "Value" = c(1,2,3,4))
test2 <- data.frame("Taxa" = c("A", "B", "C", "D"), "Value" = c(11,12,13,14))
sample_order <- data.frame("Taxa" = sample(test1$Taxa))

# Shuffled order
sample_order

# Reorder both datasets according to randomized sample order
test1
test1[order(match(test1$Taxa, sample_order$Taxa)), ]
test2
test2[order(match(test2$Taxa, sample_order$Taxa)), ]



