# ---------------------------------------
# Merritt Burch
# mbb262@cornell.edu
# 2019-07-24
# Script to visulize TASSEL GWAS results
# ---------------------------------------

# Set workdir
setwd("~/Box Sync/Cornell_PhD/labProjects/hap_gwa/July_18")

# Import data
results <- read.table("GLM_Stats_ZeaGBSv27_publicSamples_imputedV5_AGPv4-181023_Filter_Filter_MinCount2000_MinAlleleFrq01 + phenotypes_familyCovars_ClassSize10.txt",
                      header = TRUE)

# Subset down the file to be compatible with qqman
asi <- results[which(results$Trait == "ASI_BLUP_Sum0607"), c(2:4,6)]
dta <- results[which(results$Trait == "Days_To_Anthesis_BLUP_Sum0607"), c(2:4,6)]
dts <- results[which(results$Trait == "Days_To_Silk_BLUP_Sum0607"), c(2:4,6)]

# Load the library
library(qqman)

# Make the Manhattan plot on the gwasResults dataset
manhattan(asi, chr="Chr", bp="Pos", snp="Marker", p="p", annotatePval = 5.3689e-10, 
          suggestiveline = -log10(0.05), genomewideline = -log10(0.05/nrow(asi)))

manhattan(dta, chr="Chr", bp="Pos", snp="Marker", p="p", annotatePval = 2.6531e-10, 
          suggestiveline = -log10(0.05), genomewideline = -log10(0.05/nrow(dta))) 

manhattan(dts, chr="Chr", bp="Pos", snp="Marker", p="p", annotatePval = 2.6201e-10, 
          suggestiveline = -log10(0.05), genomewideline = -log10(0.05/nrow(dts)))

# ggman
library(ggman)

png("asi_nam_familyIDs.png",width = 5, height = 3.5, units = "in", res = 500)
ggman(asi, chrom="Chr", bp="Pos", snp="Marker", pvalue="p", pointSize = 0.3, lineColour = "black",
      title = "Anthesis-Silking-Interval in NAM (Family as factor)")
dev.off()

png("dta_nam_familyIDs.png",width = 5, height = 3.5, units = "in", res = 500)
ggman(dta, chrom="Chr", bp="Pos", snp="Marker", pvalue="p", pointSize = 0.3, lineColour = "black",
      title = "Days to Anthesis in NAM (Family as factor)")
dev.off()

png("dts_nam_familyIDs.png",width = 5, height = 3.5, units = "in", res = 500)
ggman(dts, chrom="Chr", bp="Pos", snp="Marker", pvalue="p", pointSize = 0.3, lineColour = "black",
      title = "Days to Silking in NAM (Family as factor)")
dev.off()


# PCA results
pcas <- read.table("MDS_eigenvalues.txt", header = TRUE)
mds <- read.table("MDS_PCA_results.txt", header = TRUE)
plot(mds$PC1,mds$PC2)
plot(mds$PC2,mds$PC3)
plot(mds$PC4,mds$PC5)
plot(mds$PC1,mds$PC5)
