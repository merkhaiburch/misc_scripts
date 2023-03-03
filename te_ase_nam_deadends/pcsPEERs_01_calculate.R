# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-01-08
# Updated... 2023-01-08
#
# Description:
# Calculate PCs and PEERs on NAM inbreds and hybrids
# ------------------------------------------------------------------------------

# Set memory
options(java.parameters = c("-Xmx30g"))
rTASSEL::startLogger(fullPath = NULL, fileName = "~/Desktop/rtassel_debug.txt")

# Load packages
library(dplyr)
library(rTASSEL)
library(peer) # Use the one installed on cbsu

# Tried installing on local computer but cmake is cursed


# ------------------------------------------------------------------------------
## Calculate PCs
# ------------------------------------------------------------------------------

# Load in genotype table
tasGenoVCF <- rTASSEL::readGenotypeTableFromPath(
  path = "~/Desktop/TutorialData/mdp_genotype.hmp.txt"
)

# Calculate PCs
pcaRes <- pca(tasGenoHMP)

# Calculate MDS
mdsRes <- mds(tasDist)


# ------------------------------------------------------------------------------
## Calculate PEERs
# ------------------------------------------------------------------------------

# PEER variables 
covariates=TRUE
PEER_number=25
tissues <- c("L3Tip", "L3Base", "FlagLeafFlowering", "GrowingPoint")

# Path to expression matrices
df_path <- "/workdir/mbb262"

# Iterate through each tissue, calculate peers
for (i in c("GRoot","GShoot","Kern","L3Base","L3Tip","LMAD","LMAN")){
  
  # Find expression matrix with matching tissue
  expression_df <- read_delim(file=paste(basedir, "some_pattern", i,".txt", sep=""), col_names=T, delim=" ")
  
  # Format df
  Taxa_list_w_tissue=expression_df[,1:10]
  colnames(expression_df)[6]="ID"
  expression_df_only_ID=expression_df[-c(1:5, 7:10)]
  Taxa_list=as.data.frame(expression_df_only_ID$ID)
  colnames(Taxa_list)="HMP32Name"
  
  # Account for the Genetic PCs when calculating the PEERs
  if (covariates==TRUE){
    covdir="/media/kak268/B_2TB_Internal/Genotypes/HM32/LDKNNi/PCs/"
    covs=read.table(paste(covdir, "MDS_5coords_PEER_fmt_merged_flt_c1-10.hmp321.onlyRNAset_MAFover005.KNNi.txt", sep=""), header=T )
    
    #repeat rows from the covariate table to match the phenotype length so that the taxa are repeated which have multiple expression phenotypes or tissues
    covs_w_rep=merge(x=Taxa_list, y=covs, by="HMP32Name")
  }
  
  # Calculate PEER factors
  model=PEER()
  
  PEER_setPhenoMean(model, as.matrix(expression_df_only_ID[,2:ncol(expression_df_only_ID)]))
  dim(PEER_getPhenoMean(model))
  
  PEER_setNk(model, PEER_number)
  
  # Adjust for covariates
  PEER_setCovariates(model, as.matrix(covs_w_rep[, 2:ncol(covs_w_rep)]))
  dim(PEER_getCovariates(model))
  
  num_of_PCs=ncol(PEER_getCovariates(model))
  
  PEER_update(model) # if you accidentally import a different number of cov's and try to run PEER then R will crash 
  PEER_plotModel(model)
  
  
  #Write the PEER factors to a file
  covs_w_PEER_factors=as.matrix(PEER_getX(model))
  
  cov_name_list="HMP32Name"
  for (x in 1:ncol(covs_w_PEER_factors)){ #add row names of cov or PEER factor
    if (x<=num_of_PCs){
      cov_name_list=c(cov_name_list, paste("geneticPC",x,sep=""))
    }
    else # subtract num_of_PCs from x below becasue there are 5 covariates ahead of the peer factors
      cov_name_list=c(cov_name_list, paste("PEER_factor",x-num_of_PCs,sep=""))
  }
  
  covs_w_PEER_factors_and_TaxaNames=cbind(Taxa_list, covs_w_PEER_factors)
  colnames(covs_w_PEER_factors_and_TaxaNames)=cov_name_list
  
  covs_w_PEER_factors_and_TaxaNames_w_TissueNames=cbind(Taxa_list_w_tissue, covs_w_PEER_factors_and_TaxaNames[,2:ncol(covs_w_PEER_factors_and_TaxaNames)])
  
  # Write results to file
  outputdir=""
  data.table::fwrite(covs_w_PEER_factors_and_TaxaNames_w_TissueNames, 
              file=paste(outputdir, i,"/", "MDSPCs_from_HMP321_w_", PEER_number,"_PEER_factors_", i,"_complete_taxa_and_tiss_names.csv", sep=""))
}
