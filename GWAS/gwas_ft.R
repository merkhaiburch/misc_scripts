# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-05-14 
# Updated... 2022-05-14

# Description 
# Check mapping of traits in the 282 for QC reasons
# ---------------------------------------------------------------

# Load in packages
library(data.table)
library(dplyr)
library(ggplot2)


# Try looking at buckler phenotypes
flowering_phenos <- data.table::fread("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/data/phenos/2021_09_final_update/goodman/02_goodman_phenotypes_Buckler_2009.txt")
colnames(flowering_phenos)[1] <- "Trait"
cor(flowering_phenos[,2], flowering_phenos[,3])

flowering_phenos_2 <- data.table::fread("~/git_projects/pleiotropy/data/all_Assoc282_Phenos.csv")
flowering_phenos_2 <- flowering_phenos_2 %>% select(Geno_Code, hapmap_names, DTA_BLUP_Buckler_2009_goodman, DTS_BLUP_Buckler_2009_goodman, ASI_BLUP_BLUP_Buckler_2009_goodman)
temp <- merge(flowering_phenos, flowering_phenos_2, by.x = "Trait", by.y = "hapmap_names")

# Load in original file
og_flowering <- data.table::fread("~/Box Sync/Cornell_PhD/phenotype_mapping_gwas/buckler_2009_science_wholeplant_floweringtime/NAMSum0607FloweringTraitBLUPsAcross8Envs.csv")
colnames(og_flowering)[1] <- "Trait"

# Merge two files together
temp2 <- merge(temp, og_flowering, by.x = "Geno_Code", by.y = "Trait")

# The phenotypes are the same, the problem isn't here!!
cor(temp2[,3], temp2[,13])
cor(temp2[,4], temp2[,14])
cor(temp2[,15], temp2[,15])


# -------------------------------------------------------------------------------------------

# Load in flowering time results from buckler 2009 NAM -----------------------
buckler_files <- list.files("~/Box Sync/Cornell_PhD/labProjects/debugging/",
                            pattern = "*Buckler_2009*", full.names = TRUE)
buckler2009 <- rbindlist(lapply(buckler_files, data.table::fread))

# Subset into dta, dts, and asi results
dta <- buckler2009 %>% filter(Trait == "DTA_BLUP_Buckler_2009_goodman")
dts <- buckler2009 %>% filter(Trait == "DTS_BLUP_Buckler_2009_goodman")
asi <- buckler2009 %>% filter(Trait == "ASI_BLUP_BLUP_Buckler_2009_goodman")

# Plot
plot_gwas_qq(dta, "DTA ~ SNP + 3 PCs (ames to 282)", "dta_goodman")
plot_gwas_qq(dts, "DTS ~ SNP + 3 PCs (ames to 282)", "dts_goodman")
plot_gwas_qq(asi, "ASI ~ SNP + 3 PCs (ames to 282)", "asi_goodman")


# Load in kernel color, shape, and flowering traits from Romay 2013 ----------
romay_files <- list.files("~/Box Sync/Cornell_PhD/labProjects/debugging/",
                            pattern = "*Romay_2013*", full.names = TRUE)
romay2013 <- rbindlist(lapply(romay_files, data.table::fread))

# Subset into sweet and flowering
sweet <- romay2013 %>% filter(Trait == "Sweet_BLUP_Romay_2013_goodman")
gdd_dts <- romay2013 %>% filter(Trait == "GDD_DTS_BLUP_Romay_2013_goodman")

# Plot
plot_gwas_qq(sweet, "sweet ~ SNP + 3 PCs (ames to 282)", "sweet_goodman")
plot_gwas_qq(gdd_dts, "GDD DTS ~ SNP + 3 PCs (ames to 282)", "gdd_dts_goodman")


# ------------------------------------------------------------------------------------
# Test mapping with PCs and MDS
# ------------------------------------------------------------------------------------

# Load in subsampled vcf file
# scp mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/mbb262/genotypes/goodman282/filtered_for_pca/goodman282_by_ames_sites_for_pca_allChroms.vcf
tasGenoVCF

# Calculate 5 PCs on SNPs
pcaRes <- pca(tasGenoVCF)
colnames(pcaRes) <- c("Taxa", "PC1", "PC2", "PC3", "PC4", "PC5")

# Calculate 5 MDS
tasDist <- distanceMatrix(tasObj = tasGenoVCF)
mdsRes <- mds(tasDist)

# Load in phenotype files --> buckler 2009 and Romay 2013

# -----------------------------------------------------------------------------------
# Map with 5 PCs
# -----------------------------------------------------------------------------------

# Setting memory
options(java.parameters = c("-Xmx250g"))
numThreads <- 55
p_value <- 0.001

# Load package
library(rTASSEL)
library(dplyr)

# Set directory
setwd("/workdir/mbb262/results/goodman_all/")

# Start logging
rTASSEL::startLogger(fullPath = NULL, fileName = NULL)

# Make directories for data and results
# mkdir -p /workdir/mbb262/phenotypes/goodman
# mkdir -p /workdir/mbb262/genotypes/goodman
# mkdir -p /workdir/mbb262/results/goodman_all/processed_results
 
# All phenotypes and covariates are in 
# scp mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/mbb262/phenotypes/0_curated_phenotypes/goodman_buckler/tasselized_for_fa/*.txt
# scp mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/mbb262/phenotypes/0_curated_phenotypes/goodman_buckler/tasselized_for_fa/* /workdir/mbb262/phenotypes/goodman
# scp mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/mbb262/genotypes/goodman282/*.gz /workdir/mbb262/genotypes/goodman

# List of phenotypes
phenos <- c("02_goodman_phenotypes_Buckler_2009.txt","09_goodman_phenotypes_Romay_2013.txt")

# Vcf path
vcfpath <- "/workdir/mbb262/genotypes/goodman/"
phenopath <- "/workdir/mbb262/phenotypes/goodman/"

# Where to output results
setwd("/workdir/mbb262/results/goodman_all") 

for (i in seq_len(10)){
  message("I am on chromosome ", i)
  
  # Load in genotype table
  vcf <-  rTASSEL::readGenotypeTableFromPath(path = paste(vcfpath,"hmp321_282_agpv4_merged_chr", i, "_imputed_goodman282.vcf.gz", sep = ""),
                                             keepDepth = FALSE)
  print(vcf)
  
  # Iterate through phenotypes
  for (j in seq_along(phenos)){
    message("I am on phenotype file ", phenos[j])
    
    # Load in phenotype file
    phenotype <- data.table::fread(paste(phenopath, phenos[j], sep = ""), nThread = numThreads)
    colnames(phenotype)[1] <- "Taxa"
    
    # Intersect phenotypes and PCs
    pheno_gpcs <- merge(x = phenotype, y = pcaRes)
    print(str(pheno_gpcs))
    
    # Tasselize merged phenotypes + global PCs + window PCs
    tasPhenoDF <- rTASSEL::readPhenotypeFromDataFrame(
      phenotypeDF = pheno_gpcs,
      taxaID = "Taxa",
      attributeTypes = c(rep("data", ncol(phenotype)-1), rep("covariate", ncol(pcaRes)-1)))
    
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
      siteMinCount = 100)
    print(tasGenoPhenoFilt)
    
    # Run fast association, write files to disk.
    rTASSEL::assocModelFitter(
      tasObj = tasGenoPhenoFilt,
      formula = . ~ .,
      fitMarkers = TRUE,
      kinship = NULL,
      fastAssociation = TRUE,
      maxP = p_value,
      maxThreads = numThreads,
      outputFile = paste("chrom_", i, "_fast_assoc_results_5pcs_", phenos[j], sep = "")
    )
  }
}






