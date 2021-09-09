# ---------------------------------------
# Merritt Burch
# mbb262@cornell.edu
# 2019-06-27
# Script to remove NAM taxa duplicates,
# take the first match only
# ---------------------------------------

# Set directory
setwd("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/data/data_summaries")

# Load packages
library(tidyverse)

# Read data and tibble it
taxa <- as_tibble(read.table("taxa_NAM.txt", header = T))
taxaSum <- as_tibble(read.table("NAM_TaxaSummary.txt", header = TRUE))
overallSum <- as_tibble(read.table("NAM_overallSummary.txt", header = TRUE))


# ---------------
# Filter by name
# ---------------

# Copy column names to work with
taxaSum2 <- taxaSum
taxaSum2$copiedTaxaNames <- taxaSum2$Taxa_Name

# Remove numbers at the end
# Next time try with [0-9]{n}
taxaSum2$copiedTaxaNames <- gsub("*:[0-9]{9}", "", taxaSum2$copiedTaxaNames)

# Get just NAM RIL names
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
taxaSum2$rilNums <- substrRight(taxaSum2$copiedTaxaNames, 4)

# Remove RILs > 200
# From 5874 to 5646
taxaSum2$rilNums <- as.numeric(taxaSum2$rilNums)
taxaSum2 <- taxaSum2[which(taxaSum2$rilNums < 201),]


# Get RIL line with the most data
aa <- taxaSum2[order(taxaSum2$copiedTaxaNames, taxaSum2$Gametes_Missing), ] #sort by id and value
taxaSum2 <- aa[ !duplicated(aa$copiedTaxaNames), ] 

# See the distribution of number of rils and ril numbers per family
tmp <- as.data.frame(as.factor(taxaSum2$rilNums))
temp <- table(tmp)

# Save files
write.table(taxaSum2$Taxa_Name, "NAM_taxa_main200_mostData", row.names = FALSE, quote = FALSE)


# -------------------------------------------
# Merge short nam names with long NAM names
# -------------------------------------------

# Import phenotypes
phenos <- read.table("NAM_phenotypes.txt", header = TRUE, fill = TRUE)

# Merge dataframes
together <- merge(phenos, taxaSum2, by.x = "Geno_Code", by.y = "copiedTaxaNames")

# Make+write New phenotype file with long NAM names
tmp <- together[,c(6,2,3,4)]
write.table(tmp, "longNAMPheno.txt", row.names = FALSE, quote = FALSE, sep = "\t")

# Make column for covariates
phenoCovar <- read.table("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/data/longNAMPheno.txt",
                         header = TRUE)

# Pull out taxa information
copiedTaxaNames <- data.frame(phenoCovar[,1])
copiedTaxaNames[,1] <- gsub("*E[0-9]{4}:[0-9]{9}", "", copiedTaxaNames[,1])
# copiedTaxaNames <- cbind(phenoCovar[,1], copiedTaxaNames)
colnames(copiedTaxaNames) <- c("<Trait>", "Q1")
old <- read.table("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/July_18/longNAMPheno_phenotype copy.txt", header = TRUE)
together <- cbind(old, copiedTaxaNames[,1])
write.table(together, "familyCovars.txt", row.names = FALSE, quote = FALSE, sep = "\t")

# -----------------
# Do same for IBMs
# -----------------

# Load packages
library(tidyverse)

# Import all IBM taxa
ibm <- as_tibble(read.delim("IBM_taxaSummary.txt"))

# Copy column names to work with
ibm$copiedTaxaNames <- ibm$Taxa.Name

# Remove numbers at the end
# Next time try with [0-9]{n}
ibm$copiedTaxaNames <- gsub("*:[0-9]{9}", "", ibm$copiedTaxaNames)

# Remove taxa with more than 30% missing data
ibm <- ibm[which(ibm$Proportion.Missing < 0.3),]

# Keep RIL lines with the most data
aa <- ibm[order(ibm$copiedTaxaNames, +ibm$Proportion.Heterozygous, ibm$Proportion.Missing), ] #sort by id and value
ibm <- aa[ !duplicated(aa$copiedTaxaNames), ] 

# Save files
write.table(ibm$Taxa.Name, "IBM_mostData.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)

